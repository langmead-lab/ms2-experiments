#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
# Fig Rep series B counting run for one new chr22 rung, as the freeze a1-chr22-ladder cells ran:
# fig_rep.py run with ms2 2885ee8 ms2-replog (one thread, interleave 32) on the simulated APR008
# chr22 PacBio reads, then fig_rep.py summarize and fig_rep_violin.py --plus-one.
# Usage: sbatch -J orbitms-figrep-run-CELL ... run_b.sbatch CELL   (CELL like chr22-h200)
set -euo pipefail
cell=$1
W=/scratch4/blangme2/ms-bakery/figrep-axes
B=/scratch16/blangme2/ms-bakery
C=$B/bench/campaign-20261012
READS=$B/results/freeze-d13/inputs/chr22-sim-pacbio-60mb-v2/reads.fa
I=$W/build/$cell
name=chr22-B-${cell#chr22-}
cd $W/runs
echo "8e5897988e4744c12a8f7fe9bf2734a38bb7ed3cdb9df40013370040c9da177f  $READS" | sha256sum -c
test -e $name/summary.json || /usr/bin/time -v python3 $C/fig_rep.py run --binary $W/tools-root/tools/ms2-2885ee8/ms2-replog \
    --reads $READS --destination $name --interleave 32 \
    --index-psi $I/index.tms.psi.orbitidx --index-phi $I/index.tms.phi.orbitidx --index-dual $I/index.tms.full.orbitidx
python3 $C/fig_rep.py summarize $name > $name.summarize.log 2>&1
PYTHONPATH=$C python3 $W/code/fig_rep_violin.py --plus-one violins-$name.json $name
echo "done $name"
