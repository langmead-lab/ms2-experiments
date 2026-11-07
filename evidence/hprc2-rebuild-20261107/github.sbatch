#!/bin/bash
#SBATCH -J orbitms-hprc2-rebuild-github
#SBATCH -p shared,parallel
#SBATCH -A blangme2
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH -o %x-%j.out
# Public-clone check: stage tools of rebuild_hprc2.sh from a fresh GitHub clone of ms2-experiments main, cloning
# TeraTools and runperm from GitHub (no TERATOOLS_REPO or RUNPERM_REPO override). The TeraLCP it builds then runs
# on the deposited mtb152 FMD, as in data/prep/mtb152_split.sh, and its outputs are compared with the deposit.
set -euo pipefail
W=/tmp/orbitms-hprc2-github-$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
df -h /tmp; hostname; gcc --version | head -1
unset TERATOOLS_REPO RUNPERM_REPO
export MS2EXP_DATA=$W/data THREADS=16
D=$MS2EXP_DATA
git clone -q https://github.com/langmead-lab/ms2-experiments.git $W/repo
cd $W/repo; git rev-parse HEAD | tee $SLURM_SUBMIT_DIR/github-repo-head.txt
python3 -m ms2exp.fetch --data $D tools/ropebwt3-cfad8475/ropebwt3 tools/runtime/libgcc_s.so.1 tools/runtime/libgomp.so.1 \
    ropebwt3/mtb152/index.fmd
chmod +x $D/tools/ropebwt3-cfad8475/ropebwt3
/usr/bin/time -v data/prep/rebuild_hprc2.sh $W/work tools 2> >(tee $SLURM_SUBMIT_DIR/github-tools.time >&2)
git -C $W/work/teratools remote get-url origin; git -C $W/work/runperm remote get-url origin
TL=$W/work/teratools/src/TeraLCP/TeraLCP; F2R=$W/work/runperm/src/ms/tools/fmd2rlbwt/fmd2rlbwt
want() { awk -F'\t' -v n="$1" '$1 == n {print $3}' data/inputs.tsv; }
report() {  # NAME FILE
    local got; got=$(sha256sum "$2" | cut -d' ' -f1)
    if [ "$got" = "$(want $1)" ]; then echo "equal   $1 $got"; else echo "DIFFER  $1 $got expected $(want $1)"; fi
}
R=$SLURM_SUBMIT_DIR/github-result.txt; : > $R
M=$W/mtb152; mkdir -p $M/tmp
$F2R $D/ropebwt3/mtb152/index.fmd $M/mtb152 > $M/fmd2rlbwt.txt
$TL -f fmd -i $D/ropebwt3/mtb152/index.fmd -t $M/tmp/t -p 16 -ominima $M/mtb152.min -ophi $M/mtb152.phi > $M/teralcp.stdout 2>&1
for f in bwt.heads bwt.len min phi; do report precursors/mtb152/$f $M/mtb152.$f | tee -a $R; done
