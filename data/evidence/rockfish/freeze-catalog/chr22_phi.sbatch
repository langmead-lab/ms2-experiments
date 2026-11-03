#!/bin/bash
#SBATCH -p shared
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=90
#SBATCH -J orbitms-freeze-chr22phi
#SBATCH -o /scratch16/blangme2/ms-bakery/freeze-catalog/chr22-phi/%x-%j.out
# TeraLCP -ophi precursors for the FMD-matched chr22 ladder cells, which the freeze build jobs need:
# their BWTs have one terminator per string, so tms-build takes phi from -ophi samples (--phi) and the
# LCP summaries from the renumbered minima (--minima). Then build phi and full with the phi-promote
# binary on chr22-h10 and check that psi, phi and dual give one output on 2,000 reads.
set -euo pipefail
TERALCP=/vast/blangme2/hprc2-ms/bin/TeraLCP
MS=/scratch16/blangme2/ms-bakery/results/campaign-20261012/bin-explore/2134460/ms-bench
OUT=/scratch16/blangme2/ms-bakery/freeze-catalog/chr22-phi
for CELL in ${CELLS:-chr22-h10 chr22-h100}; do
  SRC=/scratch16/blangme2/ms-bakery/results/ladder/$CELL/ropebwt3/index.fmd
  D=$OUT/$CELL; mkdir -p $D/tmp; cd $D
  if [ ! -s fmd.phi ]; then
    /usr/bin/time -v $TERALCP -f fmd -i $SRC -t $D/tmp/tmp -p 8 -ophi $D/fmd.phi.partial -ominima $D/teralcp.min > teralcp.stdout 2> teralcp.time
    mv fmd.phi.partial fmd.phi
  fi
  cmp teralcp.min /vast/blangme2/hprc2-ms/dev/fmdin/$CELL/fmd.min && echo "$CELL: minima match the renumbered fmdin fmd.min"
  ls -la fmd.phi; sha256sum fmd.phi
done
T=${TEST:-chr22-h10}; D=$OUT/$T; P=/vast/blangme2/hprc2-ms/dev/fmdin/$T/fmd; cd $D
for L in phi full; do
  /usr/bin/time -v $MS tms-build $P.bwt.heads $P.bwt.len $L.idx --layout $L --phi fmd.phi --minima $P.min > build-$L.out 2> build-$L.time
  tail -1 build-$L.out
done
$MS tms-build $P.bwt.heads $P.bwt.len psi.idx --layout psi > build-psi.out 2>&1
python3 - <<'PY'
import random
random.seed(1); out = open('reads.fa', 'w'); total = 0
seq = []
for line in open('/scratch16/blangme2/ms-bakery/results/ladder/chr22-h10/ropebwt3/reference.fa'):
    if line.startswith('>'):
        if len(seq) > 0: break
        continue
    seq.append(line.strip()); total += len(seq[-1])
    if total > 30_000_000: break
s = ''.join(seq).upper()
i = 0
while i < 2000:
    p = random.randrange(len(s) - 300); r = list(s[p:p + 250])
    if 'N' in r: continue
    i += 1
    for j in random.sample(range(250), 5): r[j] = random.choice('ACGT')
    out.write(f'>r{i}\n{"".join(r)}\n')
PY
$MS query psi.idx reads.fa --mode psi -o psi.ms
for x in "phi phi" "full phi" "full dual" "full psi"; do set -- $x
  $MS query $1.idx reads.fa --mode $2 -o $1-$2.ms
  cmp psi.ms $1-$2.ms && echo "$1 --mode $2 matches psi"
done
echo ALLDONE
