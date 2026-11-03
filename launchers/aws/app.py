#!/usr/bin/env python3
import json
import os

import aws_cdk as cdk

from freeze_stack.freeze_stack import FreezeBenchStack

# FREEZE_CONFIG selects the config (default config.json, a copy of config.example.json). The account and
# region come from the config when it names them, else from the credentials the CDK runs with
# (CDK_DEFAULT_ACCOUNT, CDK_DEFAULT_REGION); the S3 root (bench_s3_root) always comes from the config.
with open(os.environ.get('FREEZE_CONFIG', 'config.json')) as f:
    config = json.load(f)

app = cdk.App()
# stack_name lets a second run (the thread-scaling one) deploy beside the freeze stack.
FreezeBenchStack(app, config.get('stack_name', "Ms2FreezeBenchStack"),
                 env=cdk.Environment(account=config.get('account') or os.environ.get('CDK_DEFAULT_ACCOUNT'),
                                     region=config.get('region') or os.environ.get('CDK_DEFAULT_REGION')))
app.synth()
