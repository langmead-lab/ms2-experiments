#!/bin/bash
#SBATCH -J orbitms-d13-prep-a
#SBATCH -p shared
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH -o /scratch16/blangme2/ms-bakery/results/freeze-d13/logs/%x-%j.log
# D13 read sets that slice existing truth: 60 Mb caps, Fig 3 samples, the real short set, then locate prefixes.
set -euo pipefail
X="python3 /scratch16/blangme2/ms-bakery/bench/freeze/d13_prep.py"
O=/scratch16/blangme2/ms-bakery/results/freeze-d13/inputs
$X cap /scratch16/blangme2/ms-bakery/results/campaign-20261012/inputs/human_pacbio/explore-acgt-60mb-v1 human_pacbio-ms-60mb-acgt-v2
$X cap /scratch16/blangme2/ms-bakery/results/campaign-20261012/inputs/human_ont_modern/explore-acgt-60mb-v1 human_ont_modern-ms-60mb-acgt-v2
$X cap /scratch16/blangme2/ms-bakery/results/campaign-20261012/inputs/fig_ops2/human472-apr008-pacbio-60mb-v1 human472-apr008-pacbio-60mb-v2
$X cap /scratch16/blangme2/ms-bakery/results/campaign-20261012/inputs/fig_ops2/mtb152-pacbio-60mb-v1 mtb152-pacbio-60mb-v2
$X sample /scratch16/blangme2/ms-bakery/results/campaign-20261012/inputs/fig_ops/human472-apr008-hap1-p100k-v1 human472-apr008-hap1-p100k-s600-v1 --seed 2026103001
$X sample /scratch16/blangme2/ms-bakery/results/campaign-20261012/inputs/fig_ops/mtb152-heldout16-p100k-v1 mtb152-heldout16-p100k-s600-v1 --seed 2026103002
$X short-60mb
$X cap $O/human_short-ms-60mb-acgt-v2 human_short-locate-6mb-v2 --limit 6000000
$X cap $O/human_pacbio-ms-60mb-acgt-v2 human_pacbio-locate-20mb-v2 --limit 20000000
$X cap $O/human_ont_modern-ms-60mb-acgt-v2 human_ont_modern-locate-10mb-v2 --limit 10000000
