#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
# One freeze job on one Rockfish node (BENCH_FREEZE_DESIGN.md D11, D12). Resources, partition,
# account, time and --exclusive come from the sbatch command line; submit.py chooses them per job.
# Usage: sbatch [resources] -J orbitms-freeze-JOB run_job.sbatch MANIFEST JOB_ID INPUTS_ROOT RESULTS_ROOT
# Inputs are staged to node-local /tmp, verified by sha256; there is no root, so each run evicts
# its inputs with posix_fadvise instead of dropping the page cache, and the timings are for
# planning only. Missing oracle traces are recorded (D12). A job that hits its time limit is
# resumed by submitting it again: cells with done markers are skipped.
set -uo pipefail
manifest=$(readlink -f "$1"); job=$2; inputs=$3; results=$4
# FREEZE_HOME is the ms2-experiments checkout (default: the one this script is in when submitted from it).
here=${FREEZE_HOME:-$SLURM_SUBMIT_DIR}
python=${FREEZE_PYTHON:-python3}
stage=/tmp/${USER}-freeze-${SLURM_JOB_ID:-$$}
mkdir -p "$stage"
cleanup() { rm -rf "$stage"; }
trap cleanup EXIT
df -h /tmp | tail -1
numactl --hardware | head -3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PYTHONPATH="$here${PYTHONPATH:+:$PYTHONPATH}" "$python" -m ms2exp.run.runner "$job" --manifest "$manifest" --inputs-root "$inputs" --results-root "$results" \
    --stage "$stage" --platform rockfish --record-traces &
pid=$!
# Forward Slurm's termination signal so the runner's current command dies with the job.
trap 'kill -TERM "$pid" 2>/dev/null' TERM INT
wait "$pid"; status=$?
while kill -0 "$pid" 2>/dev/null; do wait "$pid"; status=$?; done
exit $status
