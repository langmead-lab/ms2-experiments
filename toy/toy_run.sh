#!/bin/bash
# Run the Rockfish jobs of the toy manifest on this node (devlangmead1), then collect.
# Usage: toy_run.sh TOY_DIR [JOB_ID ...]   (after toy_prep.sh TOY_DIR); with no job ids, every job runs.
set -uo pipefail
T=$(readlink -f "$1")
shift
REPO=$(dirname "$(dirname "$(readlink -f "$0")")")
PY="env PYTHONPATH=$REPO${PYTHONPATH:+:$PYTHONPATH} python3"
M=$T/manifest-toy.json
stage=/tmp/${USER}-ms2exp-toy
status=0
jobs=${*:-$(python3 -c "import json; print(' '.join(j['id'] for j in json.load(open('$M'))['jobs']))")}
for job in $jobs; do
    $PY -m ms2exp.run.runner $job --manifest $M --inputs-root $T/inputs-root --results-root $T/results \
        --stage $stage --platform rockfish --record-traces || status=1
done
$PY -m ms2exp.collect --manifest $M --results-root $T/results --out $T/collected --platform rockfish || status=1
echo "toy run status $status"
exit $status
