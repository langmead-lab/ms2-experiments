#!/bin/bash
# Publish a freeze run to S3 for the AWS boxes: the code tarball, the manifest, and every input.
# Run on Rockfish, where the inputs tree is, after the Rockfish dry run passed (D12).
# Usage: publish.sh MANIFEST INPUTS_ROOT S3_ROOT [TRACES_DIR] [--dry-run]
#   S3_ROOT is config.json's bench_s3_root; inputs go to S3_ROOT/inputs/<logical name>,
#   TRACES_DIR is RESULTS_ROOT/MANIFEST_ID/traces from the Rockfish runs. With INPUTS_S3 set (the
#   config's inputs_s3_root), inputs go there instead; ones already present are skipped.
# Needs credentials for the fleet account's freeze bucket, for example AWS_PROFILE=index-zone-ec2
# after `aws sso login --profile index-zone-ec2 --use-device-code`. This uploads; it is never run
# by the workflow itself.
set -euo pipefail
manifest=$(readlink -f "$1"); inputs=$2; root=${3%/}; traces=${4:-}; dry=${5:-}
repo=$(cd "$(dirname "$0")/../.." && pwd)
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
# The boxes run the ms2exp package and launchers/aws/worker.py from the repository root.
(cd "$repo" && tar czf "$tmp/freeze-code.tar.gz" --exclude=__pycache__ --exclude=.venv --exclude=cdk.out \
    ms2exp launchers/__init__.py launchers/aws/__init__.py launchers/aws/worker.py)
run() { if [ "$dry" = --dry-run ]; then echo "$@"; else "$@"; fi; }
run aws s3 cp --only-show-errors "$tmp/freeze-code.tar.gz" "$root/code/freeze-code.tar.gz"
run aws s3 cp --only-show-errors "$manifest" "$root/$(python3 -c "import json;print(json.load(open('$(dirname "$0")/${FREEZE_CONFIG:-config.json}'))['manifest'])")"
PYTHONPATH="$repo" python3 -m ms2exp.catalog publish --manifest "$manifest" --root "$inputs" --s3 "${INPUTS_S3:-$root/inputs}" \
    ${traces:+--traces "$traces"} ${dry:+--dry-run}
