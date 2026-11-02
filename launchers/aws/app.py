#!/usr/bin/env python3
import json
import os

import aws_cdk as cdk

from freeze_stack.freeze_stack import FreezeBenchStack

# FREEZE_CONFIG selects another config, such as config-smoke.json for a one-box toy run. The config names
# the AWS account and the S3 root (bench_s3_root), so neither is written into the code.
with open(os.environ.get('FREEZE_CONFIG', 'config.json')) as f:
    config = json.load(f)

app = cdk.App()
# stack_name lets a second run (the thread-scaling one) deploy beside the freeze stack.
FreezeBenchStack(app, config.get('stack_name', "Ms2FreezeBenchStack"),
                 env=cdk.Environment(account=config['account'], region=config['region']))
app.synth()
