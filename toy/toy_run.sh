#!/bin/bash
# Run every job of the toy manifest directly on this machine (no Slurm), then collect.
# Usage: toy_run.sh TOY_DIR [JOB_ID ...]   (after toy_prep.sh TOY_DIR); with no job ids, every job runs.
# STAGE is the local disk the runner copies inputs to, about 6 GB. It defaults to TOY_DIR/stage,
# which is removed once every job has passed unless KEEP_STAGE=1; a STAGE given is kept.
# Without numactl the runner runs unbound (--no-numa). When transparent huge pages are not set to
# always, as the protocol times runs, the runner runs anyway (--any-thp) and records the override;
# the toy's timings are for orientation only.
# The platform is called rockfish for any host that is not an AWS box of the launcher: it selects
# the cells that run off AWS, such as the counting cells.
set -uo pipefail
T=$(readlink -f "$1")
shift
REPO=$(dirname "$(dirname "$(readlink -f "$0")")")
PY="env PYTHONPATH=$REPO${PYTHONPATH:+:$PYTHONPATH} python3"
M=$T/manifest-toy.json
stage=${STAGE:-$T/stage}
flags=()
if ! command -v numactl > /dev/null; then
    echo "numactl not found: running without NUMA binding" >&2
    flags+=(--no-numa)
fi
thp=$(sed 's/.*\[\(.*\)\].*/\1/' /sys/kernel/mm/transparent_hugepage/enabled 2>/dev/null || echo unknown)
if [ "$thp" != always ]; then
    echo "transparent huge pages are $thp, not always: running anyway (--any-thp)" >&2
    flags+=(--any-thp)
fi
status=0
jobs=${*:-$(python3 -c "import json; print(' '.join(j['id'] for j in json.load(open('$M'))['jobs']))")}
for job in $jobs; do
    $PY -m ms2exp.run.runner $job --manifest $M --inputs-root $T/inputs-root --results-root $T/results \
        --stage $stage --platform rockfish --record-traces ${flags[@]+"${flags[@]}"} || status=1
done
$PY -m ms2exp.collect --manifest $M --results-root $T/results --out $T/collected --platform rockfish || status=1
if [ $status = 0 ] && [ -z "${STAGE:-}" ] && [ "${KEEP_STAGE:-0}" != 1 ]; then
    rm -rf -- "$T/stage"
fi
echo "toy run status $status"
exit $status
