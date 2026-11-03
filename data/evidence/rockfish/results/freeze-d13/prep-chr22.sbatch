#!/bin/bash
#SBATCH -J orbitms-d13-chr22sim
#SBATCH -p shared
#SBATCH -c 16
#SBATCH --mem=24G
#SBATCH --time=6:00:00
#SBATCH -o /scratch16/blangme2/ms-bakery/results/freeze-d13/logs/%x-%j.log
# D13.1: APR008 hap1 chr22 (q-arm scaffold), then simulated short, PacBio and ONT reads with truth on chr22-h10 and chr22-h100.
set -euo pipefail
test -e /scratch16/blangme2/ms-bakery/results/freeze-d13/inputs/haplotype-APR008-hap1-chr22-v1 || python3 /scratch16/blangme2/ms-bakery/bench/freeze/d13_prep.py apr008-chr22 --scaffold JBHMUJ010000002.1 --evidence "minimap2 2.28 -x asm5 of every APR008 hap1 scaffold against CHM13 v2.0 chr22 (job 32742987, freeze-d13/apr008-chr22): JBHMUJ010000002.1 covers chr22 17.4-51.3 Mb (33.85 Mb in alignments of at least 500 kb, 99.5 percent of the scaffold) and no other scaffold aligns there; chr22 0-17.4 Mb (short arm and centromere) is spread over at least 10 scaffolds shared with other acrocentrics and is not included"
for s in short pacbio ont; do python3 /scratch16/blangme2/ms-bakery/bench/freeze/d13_prep.py simulate --set $s --threads 16; done
