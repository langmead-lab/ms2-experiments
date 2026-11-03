#!/bin/bash
#SBATCH -p bigmem
#SBATCH -A blangme2_bigmem
#SBATCH --qos=qos_bigmem
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=120
#SBATCH -J orbitms-fmdin-A
#SBATCH -o /vast/blangme2/hprc2-ms/dev/fmdin/logs/%x-%j.out
# Stage A for one ladder cell: ropebwt3 stat, FMD runs dump, TeraLCP -f fmd minima, renumbering.
set -euo pipefail
CELL=${CELL:-chr22-h10}
SRC=/scratch16/blangme2/ms-bakery/results/ladder/$CELL/ropebwt3/index.fmd
W=/vast/blangme2/hprc2-ms/dev/fmdin
OUT=$W/$CELL
REP=/vast/blangme2/hprc2-ms/reports/fmdin
MON="python3 /vast/blangme2/hprc2-ms/scripts/monitor.py"
TOOL=$W/orbit/src/ms/tools/fmd2rlbwt/fmd2rlbwt
TERALCP=$W/TeraTools-bce9ea7/src/TeraLCP/TeraLCP
RB3=/scratch16/blangme2/ms-bakery/ropebwt3/ropebwt3
mkdir -p $OUT $REP
cd $OUT
md5sum $SRC > fmd.md5
$MON --label $CELL-stat --out $REP/$CELL-stat.monitor.tsv --watch $OUT -- $RB3 stat $SRC > stat.txt
$MON --label $CELL-dump --out $REP/$CELL-dump.monitor.tsv --watch $OUT -- $TOOL $SRC $OUT/fmd > dump.txt
mkdir -p $OUT/tl-tmp
$MON --label $CELL-teralcp-minima --out $REP/$CELL-teralcp-minima.monitor.tsv --watch $OUT -- \
  $TERALCP -f fmd -i $SRC -t $OUT/tl-tmp/tmp -p ${SLURM_CPUS_PER_TASK:-8} -ominima $OUT/teralcp.min > teralcp.stdout 2> teralcp.stderr
$MON --label $CELL-renumber-minima --out $REP/$CELL-renumber-minima.monitor.tsv --watch $OUT -- \
  $TOOL --minima-in $OUT/teralcp.min --minima-out $OUT/fmd.min $SRC $OUT/fmd > renumber.txt 2> renumber.stderr
ls -la $OUT
