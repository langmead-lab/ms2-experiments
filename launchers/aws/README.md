# AWS launcher

This directory runs a benchmark manifest on AWS spot instances ("boxes"). It is an AWS CDK app (CDK is AWS's tool for defining cloud resources in code) that keeps a group of metal instances running. Each box claims one job at a time from the manifest, runs it with `ms2exp.run.runner`, writes the results to S3, and shuts itself down when no work is left. Both freeze manifests and the thread-scaling manifest of the paper ran this way. `../../RERUN_FREEZE.md` puts this step in context, from inputs to figures.

| File | Role |
|---|---|
| `app.py`, `cdk.json`, `requirements.txt` | The CDK app and its Python dependencies |
| `freeze_stack/freeze_stack.py` | The stack: a VPC with public subnets only, an Auto Scaling group of spot boxes, and a policy on your instance role that lets boxes use the bucket and scale the group |
| `freeze_stack/userdata.sh` | What each box runs at boot: an NVMe RAID0 at `/work`, transparent huge pages set to `always`, the code from S3, then `worker.py` |
| `worker.py` | The box loop: claim a job through an S3 conditional write, run it, keep a heartbeat, resume after a spot interruption, leave when the queue is empty |
| `publish.sh` | Uploads the code tarball, the manifest and every input to the S3 root |
| `capacity_check.py` | Read-only check of spot prices, spot placement scores and your spot vCPU quota |
| `config.example.json` | The settings; copy it to `config.json` |

## What you need

- An AWS account with a spot vCPU quota for standard instances of at least 128 per box. `capacity_check.py` prints your quota and what `box_count` needs.
- A private S3 bucket in the same region. The paper's runs stored about 621 GiB of inputs and 1.6 TiB of built indexes there.
- An existing IAM role that EC2 instances can assume. The stack adds the S3 and Auto Scaling permissions the boxes need and removes them on `cdk destroy`.
- On your own machine: the AWS CLI with credentials for the account, Node.js with the CDK CLI (`npm install -g aws-cdk`), and Python 3.10 or later.

## Settings (`config.json`)

| Key | Meaning |
|---|---|
| `region` | AWS region of the boxes and the bucket |
| `instance_types` | Instance types tried in order; every box of the paper's runs was `i4i.metal` (1 TiB of memory, 2 sockets of 32 cores, 8 NVMe disks) |
| `box_count`, `initial_boxes` | Most boxes running at once, and how many start at deploy; boxes add more as jobs become ready |
| `max_idle_minutes` | How long a box with nothing to claim waits before it leaves |
| `spot_max_price` | Highest spot price per hour, or `null` for the on-demand price |
| `max_hours_per_box` | A box powers off after this many hours whatever happens; its job resumes on a replacement |
| `root_volume_gb` | Root disk size; inputs and results go to the NVMe RAID0, not here |
| `instance_role_name` | The existing EC2 instance role |
| `bench_s3_root` | S3 prefix holding `code/`, the manifest, `inputs/` and `results/` |
| `manifest` | File name of the manifest under `bench_s3_root` |

Optional keys: `ssh_cidr` (an address range such as `203.0.113.7/32`; the stack then creates a key pair and opens SSH from it), `inputs_s3_root` (read inputs that an earlier run published) and `stack_name` (deploy a second stack beside the first).

## Running a manifest

From this directory, with `MANIFEST` the manifest file and `DATA` the directory holding every input under its logical name (`MS2EXP_DATA` after `python3 -m ms2exp.fetch`):

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt && source .venv/bin/activate
cp config.example.json config.json            # then edit region, role, bucket and manifest name
python3 capacity_check.py --boxes 5           # spot prices, placement scores, quota; changes nothing
./publish.sh MANIFEST DATA s3://YOUR-BUCKET/freeze --dry-run
./publish.sh MANIFEST DATA s3://YOUR-BUCKET/freeze
cdk bootstrap                                 # once per account and region
cdk deploy
```

`publish.sh` checks each input's sha256 before uploading it and skips objects already in S3 with the same size. Build jobs run first and publish the indexes they make; timed jobs are offered as soon as the indexes they read exist.

To follow progress, count the finished jobs (a freeze manifest has 83) and read the box logs:

```bash
ID=$(python3 -c "import json; print(json.load(open('MANIFEST'))['id'])")
aws s3 ls --recursive s3://YOUR-BUCKET/freeze/results/$ID/ | grep -c '/job.done$'
aws s3 ls s3://YOUR-BUCKET/freeze/results/$ID/boxes/
```

When every job has its `job.done`, the boxes leave on their own. Then remove the stack; this leaves the bucket and its contents:

```bash
cdk destroy
```

## Cost

TBD. The paper's runs were not costed afterwards. Before the runs, the design record estimated about 110 to 125 instance-hours for the freeze (`docs/BENCH_FREEZE_DESIGN.md`, "First-pass AWS cost estimate"); multiply by the current spot price that `capacity_check.py` prints. Inputs and indexes also cost S3 storage for as long as they stay in the bucket.
