#!/bin/bash
# Boot script of one freeze benchmark box; the stack fills in the @NAME@ placeholders.
set -x
exec > /var/log/freeze-userdata.log 2>&1
# Backstop: power off (which terminates the box) after the box's hour limit whatever happens.
shutdown -h +$((@MAX_HOURS@ * 60))
# Huge pages for every process (the protocol's transparent_hugepages); Amazon Linux defaults to madvise,
# and ms2 does not ask for them. Timed jobs check this setting before they start.
echo always > /sys/kernel/mm/transparent_hugepage/enabled
dnf install -y numactl time tar gzip unzip python3 mdadm xfsprogs nvme-cli zlib
# Build jobs of an all-AWS run compile ms2 (ms2exp.build.tools links libstdc++ statically); A1 cells need numpy (fig_rep.py).
dnf install -y git make gcc-c++ libstdc++-static zlib-devel python3-numpy
command -v aws || (cd /tmp && curl -sSL https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip \
    -o awscliv2.zip && unzip -q awscliv2.zip && ./aws/install)

# NVMe instance-store RAID0 at /work (r6idn.metal has 2 disks, i4i.metal 8).
disks=$(nvme list | awk '/Amazon EC2 NVMe Instance Storage/ {print $1}')
n=$(echo $disks | wc -w)
if [ "$n" -gt 1 ]; then
    mdadm --create /dev/md0 --run --level=0 --raid-devices=$n $disks
    dev=/dev/md0
else
    dev=$disks
fi
mkfs.xfs -f -K $dev
mkdir -p /work
mount -o noatime $dev /work

# S3 access uses the instance role's own credentials; the stack grants it the freeze bucket.
export AWS_DEFAULT_REGION=@REGION@
ROOT=@BENCH_S3_ROOT@
W=/work

# The code (the ms2exp package and the launchers) and the manifest, published by
# launchers/aws/publish.sh. Tool binaries are manifest inputs: the runner stages and verifies them
# like any other input (D11).
mkdir -p $W/code
aws s3 cp --only-show-errors $ROOT/code/freeze-code.tar.gz $W/freeze-code.tar.gz
tar xzf $W/freeze-code.tar.gz -C $W/code
aws s3 cp --only-show-errors $ROOT/@MANIFEST@ $W/manifest.json
cd $W/code
setsid env PYTHONPATH=$W/code python3 -m launchers.aws.worker --manifest $W/manifest.json --root $ROOT --inputs-root @INPUTS_S3_ROOT@ --work $W \
    --terminate-when-empty --max-idle-minutes @MAX_IDLE_MINUTES@ > $W/worker.out 2>&1 < /dev/null &
echo USERDATA_DONE
