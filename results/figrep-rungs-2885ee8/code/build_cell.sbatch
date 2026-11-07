#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
# Fig Rep series B rungs above 100 haplotypes, by the chr22-h10/h100 route of the freeze:
# ladder reference (prepare_more.py) -> ropebwt3 build -d (the ladder's ropebwt3 binary)
# -> TeraLCP 414f809 -f fmd -ophi -ominima -> fmd2rlbwt (runs) -> ms2 2885ee8 build psi, phi and full with --phi and --minima, exactly as
# freeze/make_manifest.py index_steps does for FMD_REFS.
# Usage: sbatch -J orbitms-figrep-CELL ... build_cell.sbatch CELL
set -euo pipefail
cell=$1
W=/scratch4/blangme2/ms-bakery/figrep-axes
B=/scratch16/blangme2/ms-bakery
ROPE=$B/ropebwt3/ropebwt3
TERALCP=/vast/blangme2/hprc2-ms/bin/TeraLCP
FMD2RLBWT=/vast/blangme2/hprc2-ms/bin/fmd2rlbwt
MS=$W/tools-root/tools/ms2-2885ee8/ms2-bench
REF=$W/ladder/$cell/reference.fa
OUT=$W/build/$cell
mkdir -p "$OUT"; cd "$OUT"
work=/tmp/blangme2-figrep-$SLURM_JOB_ID; mkdir -p "$work"; trap 'rm -rf "$work"' EXIT
df -h /tmp | tail -1
step() {  # label argv...
    local label=$1; shift
    test -e "$label.done" && return 0
    printf '%q ' "$@" > "$label.command"; echo >> "$label.command"
    /usr/bin/time -f '%e %U %S %M' -o "$label.time" "$@" > "$label.stdout" 2> "$label.stderr"
    touch "$label.done"
    echo "$label: $(cat "$label.time") (wall user sys max_rss_kib) $(date)"
}
step prepare python3 $W/prepare_more.py $W/plan-more.json $cell --jobs 8
step ropebwt3 $ROPE build -t ${SLURM_CPUS_PER_TASK:-8} -d -o index.fmd "$REF"
step teralcp $TERALCP -f fmd -i index.fmd -t "$work/tmp" -p ${SLURM_CPUS_PER_TASK:-8} -ophi fmd.phi -ominima teralcp.min
# TeraLCP 414f809 numbers minima by the maximal runs already (for chr22-h100 its -ominima output
# equals fmdin's renumbered fmd.min), so fmd2rlbwt writes only the runs.
rm -f fmd2rlbwt.done fmd.min
step fmd2rlbwt-runs $FMD2RLBWT index.fmd fmd
ln -sfn teralcp.min fmd.min
test -x "$MS"
step build-psi $MS build fmd.bwt.heads fmd.bwt.len index.tms.psi.orbitidx --layout psi
step build-phi $MS build fmd.bwt.heads fmd.bwt.len index.tms.phi.orbitidx --layout phi --phi fmd.phi --minima fmd.min
step build-full $MS build fmd.bwt.heads fmd.bwt.len index.tms.full.orbitidx --layout full --phi fmd.phi --minima fmd.min
sha256sum "$MS" index.fmd index.tms.*.orbitidx > sha256.txt
ls -la
echo "built $cell"
