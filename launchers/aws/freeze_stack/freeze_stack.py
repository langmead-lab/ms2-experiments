"""
Spot fleet for the ms2 freeze benchmark's timed runs (BENCH_FREEZE_DESIGN.md D4, D5, D10, D11).

An Auto Scaling group keeps box_count metal spot boxes running. Its mixed-instances policy
tries r6idn.metal first and falls back to i4i.metal (same CPU, memory and topology) when
r6idn.metal capacity is short, in any us-east-1 availability zone. A box reclaimed by spot
is replaced by the group, and the replacement resumes the interrupted job from its per-cell
done markers in S3, so the group is the orchestrator that relaunches interrupted jobs.

Each box builds an NVMe RAID0 at /work from its instance-store disks, fetches the freeze
code and the pinned tool sources from bench_s3_root, builds the tools, then runs worker.py:
claim one job from the S3 queue, run it with run_job.py (one job per box, socket 0 only),
and repeat. When the queue has nothing left to claim, the box removes itself from the group
with a decrement of the desired capacity, so the group shrinks to zero instead of replacing
it. `cdk destroy` is still the backstop.

There is no perf or likwid: the paper reports wall-clock query time, not counters.

Credentials: the standing instance role is referenced, as in the other bench stacks. The stack
attaches one policy to it that lets a box scale and leave this group and read and write the
freeze bucket, a private bucket in this account named by bench_s3_root; destroy removes that
policy and leaves the bucket.
"""

import json
import os
from pathlib import Path

from aws_cdk import (
    Stack,
    Fn,
    aws_autoscaling as autoscaling,
    aws_ec2 as ec2,
    aws_iam as iam,
    CfnOutput,
    CfnTag,
)
from constructs import Construct

HERE = Path(__file__).resolve().parent


class FreezeBenchStack(Stack):
    def __init__(self, scope: Construct, construct_id: str, **kwargs) -> None:
        super().__init__(scope, construct_id, **kwargs)

        with open(os.environ.get('FREEZE_CONFIG', 'config.json')) as f:
            config = json.load(f)

        vpc = ec2.Vpc(
            self, "FreezeVpc",
            max_azs=6,
            nat_gateways=0,
            subnet_configuration=[ec2.SubnetConfiguration(
                name="freeze-public", subnet_type=ec2.SubnetType.PUBLIC, cidr_mask=24)],
        )
        sec_group = ec2.SecurityGroup(self, "FreezeSecurityGroup", vpc=vpc, allow_all_outbound=True)
        sec_group.add_ingress_rule(ec2.Peer.ipv4(config['ssh_cidr']), ec2.Port.tcp(22), "SSH")

        key_pair = ec2.CfnKeyPair(
            self, "FreezeKeyPair", key_name="ms2-freeze-bench",
            tags=[CfnTag(key="application", value="ms2-freeze-bench")])

        instance_role = iam.Role.from_role_name(
            self, "FreezeInstanceRole", config['instance_role_name'], mutable=True)

        user_data = ec2.UserData.for_linux()
        template = (HERE/'userdata.sh').read_text()
        for key, value in {
            '@REGION@': config['region'],
            '@BENCH_S3_ROOT@': config['bench_s3_root'].rstrip('/'),
            '@MANIFEST@': config['manifest'],
            '@MAX_HOURS@': str(int(config['max_hours_per_box'])),
            '@MAX_IDLE_MINUTES@': str(config.get('max_idle_minutes', 20)),
        }.items():
            template = template.replace(key, value)
        user_data.add_commands(*template.split('\n')[1:])

        launch_template = ec2.LaunchTemplate(
            self, "FreezeLaunchTemplate",
            machine_image=ec2.MachineImage.latest_amazon_linux2023(),
            security_group=sec_group,
            key_name=key_pair.key_name,
            role=instance_role,
            user_data=user_data,
            associate_public_ip_address=True,
            instance_initiated_shutdown_behavior=ec2.InstanceInitiatedShutdownBehavior.TERMINATE,
            block_devices=[ec2.BlockDevice(
                device_name="/dev/xvda",
                volume=ec2.BlockDeviceVolume.ebs(
                    config['root_volume_gb'], volume_type=ec2.EbsDeviceVolumeType.GP3,
                    delete_on_termination=True))],
        )

        types = config['instance_types']
        group = autoscaling.AutoScalingGroup(
            self, "FreezeGroup",
            vpc=vpc,
            vpc_subnets=ec2.SubnetSelection(subnet_type=ec2.SubnetType.PUBLIC),
            min_capacity=0,
            max_capacity=config['box_count'],
            # Boxes raise the desired capacity themselves as jobs become ready (worker.py).
            desired_capacity=config.get('initial_boxes', config['box_count']),
            mixed_instances_policy=autoscaling.MixedInstancesPolicy(
                launch_template=launch_template,
                launch_template_overrides=[
                    autoscaling.LaunchTemplateOverrides(instance_type=ec2.InstanceType(t)) for t in types],
                instances_distribution=autoscaling.InstancesDistribution(
                    on_demand_base_capacity=0,
                    on_demand_percentage_above_base_capacity=0,
                    # Honors the order of instance_types: r6idn.metal first, i4i.metal as fallback.
                    spot_allocation_strategy=autoscaling.SpotAllocationStrategy.CAPACITY_OPTIMIZED_PRIORITIZED,
                    spot_max_price=(str(config['spot_max_price']) if config.get('spot_max_price') else None),
                ),
            ),
            capacity_rebalance=False,
        )
        cfn_group = group.node.default_child
        cfn_group.add_property_override('Tags', [
            dict(Key='application', Value='ms2-freeze-bench', PropagateAtLaunch=True)])

        # A box removes itself from the group once the queue is empty, and raises the group's desired
        # capacity when more jobs are ready than boxes run them.
        instance_role.add_to_principal_policy(iam.PolicyStatement(
            actions=['autoscaling:DescribeAutoScalingInstances', 'autoscaling:DescribeAutoScalingGroups'],
            resources=['*'],
        ))
        # A policy of its own, since it names the group: the launch template waits for the role's
        # default policy, and the group waits for the launch template.
        instance_role.attach_inline_policy(iam.Policy(self, "FreezeScalePolicy", statements=[iam.PolicyStatement(
            actions=['autoscaling:TerminateInstanceInAutoScalingGroup', 'autoscaling:SetDesiredCapacity'],
            resources=[Fn.sub(
                'arn:${AWS::Partition}:autoscaling:${AWS::Region}:${AWS::AccountId}:autoScalingGroup:*:'
                'autoScalingGroupName/${Name}', {'Name': group.auto_scaling_group_name})],
        )]))

        bucket, _, prefix = config['bench_s3_root'][len('s3://'):].rstrip('/').partition('/')
        instance_role.add_to_principal_policy(iam.PolicyStatement(
            actions=['s3:ListBucket'],
            resources=[f'arn:aws:s3:::{bucket}'],
        ))
        instance_role.add_to_principal_policy(iam.PolicyStatement(
            actions=['s3:GetObject', 's3:PutObject', 's3:DeleteObject'],
            resources=[f'arn:aws:s3:::{bucket}/{prefix}/*' if prefix else f'arn:aws:s3:::{bucket}/*'],
        ))

        CfnOutput(self, "KeyPairId", value=key_pair.attr_key_pair_id)
        CfnOutput(self, "GroupName", value=group.auto_scaling_group_name)
