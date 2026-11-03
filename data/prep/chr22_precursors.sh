#!/bin/bash
# RLBWT precursors of the chr22-h10 and chr22-h100 references (precursors/chr22-hN/{bwt.heads,bwt.len,min,phi}),
# from the ropebwt3 FMD of data/prep/legacy_indexes.py. The frozen files came from two Rockfish jobs:
# /vast/blangme2/hprc2-ms/dev/fmdin/stageA.sh on Oct 1 2026 (copy in data/evidence/fmdin/stageA.sh) for the
# run heads, lengths and minima, and freeze-catalog/chr22_phi.sbatch on Oct 30 (copy in
# data/evidence/rockfish/freeze-catalog/) for the phi samples. In order:
#   fmd2rlbwt FMD OUT/fmd               run heads and lengths (fmd.bwt.heads, fmd.bwt.len)
#   TERALCP_MINIMA -f fmd -ominima      LCP minima from TeraTools bce9ea7, numbered by TeraLCP's runs
#   fmd2rlbwt --minima-in/--minima-out  the same minima renumbered to fmd2rlbwt's runs (fmd.min)
#   TERALCP_PHI -f fmd -ophi            phi samples (fmd.phi); its minima must equal fmd.min
# The manifest names fmd.bwt.heads, fmd.bwt.len, fmd.min and fmd.phi bwt.heads, bwt.len, min and phi.
# Usage: data/prep/chr22_precursors.sh FMD OUTDIR
#   environment: FMD2RLBWT, TERALCP_MINIMA (TeraLCP built from TeraTools bce9ea7), TERALCP_PHI (TeraLCP
#   from TeraTools 414f809, the build of the HPRC2 and mtb152 precursors), ROPEBWT3, THREADS
set -euo pipefail
SRC=$1; OUT=$2; P=${THREADS:-8}
mkdir -p $OUT/tl-tmp $OUT/phi-tmp
cd $OUT
md5sum $SRC > fmd.md5
${ROPEBWT3:-ropebwt3} stat $SRC > stat.txt
$FMD2RLBWT $SRC $OUT/fmd > dump.txt
$TERALCP_MINIMA -f fmd -i $SRC -t $OUT/tl-tmp/tmp -p $P -ominima $OUT/teralcp.min > teralcp.stdout 2> teralcp.stderr
$FMD2RLBWT --minima-in $OUT/teralcp.min --minima-out $OUT/fmd.min $SRC $OUT/fmd > renumber.txt 2> renumber.stderr
$TERALCP_PHI -f fmd -i $SRC -t $OUT/phi-tmp/tmp -p $P -ophi $OUT/fmd.phi.partial -ominima $OUT/phi-teralcp.min \
    > phi-teralcp.stdout 2> phi-teralcp.stderr
mv fmd.phi.partial fmd.phi
cmp phi-teralcp.min fmd.min && echo "minima of the phi run equal fmd.min"
sha256sum fmd.bwt.heads fmd.bwt.len fmd.min fmd.phi
