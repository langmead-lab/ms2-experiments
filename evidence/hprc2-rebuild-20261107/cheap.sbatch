#!/bin/bash
#SBATCH -J orbitms-hprc2-rebuild-cheap
#SBATCH -p shared,parallel
#SBATCH -A blangme2
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH -o %x-%j.out
# Cheap tier: public human472.fmr.gz -> index.fmd -> bwt.heads, bwt.len, each checked against data/inputs.tsv.
# Everything large stays in a per-job directory on node-local /tmp. Usage: sbatch cheap.sbatch REPO [TERATOOLS_REPO]
set -euo pipefail
REPO=$1; export TERATOOLS_REPO=${2:-https://github.com/BenLangmead/TeraTools.git}
W=/tmp/orbitms-hprc2-cheap-$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
df -h /tmp; hostname; gcc --version | head -1
export MS2EXP_DATA=$W/data THREADS=8
cd $REPO
python3 -m ms2exp.fetch --data $MS2EXP_DATA tools/ropebwt3-cfad8475/ropebwt3 tools/runtime/libgcc_s.so.1 tools/runtime/libgomp.so.1
chmod +x $MS2EXP_DATA/tools/ropebwt3-cfad8475/ropebwt3
for s in tools fmd runs; do
  /usr/bin/time -v data/prep/rebuild_hprc2.sh $W/work $s 2> >(tee $SLURM_SUBMIT_DIR/cheap-$s.time >&2)
done
cp $W/work/fmd2rlbwt.txt $SLURM_SUBMIT_DIR/
