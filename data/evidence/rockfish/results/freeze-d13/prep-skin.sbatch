#!/bin/bash
#SBATCH -J orbitms-d13-skin
#SBATCH -p shared
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH -o /scratch16/blangme2/ms-bakery/results/freeze-d13/logs/%x-%j.log
# D13 skin set: filter first, cap at 60,000,000 bases, new truth on the human472 FMD; then its 0.5 Mb locate prefix.
set -euo pipefail
X="python3 /scratch16/blangme2/ms-bakery/bench/freeze/d13_prep.py"
$X skin-60mb --threads 16
$X cap /scratch16/blangme2/ms-bakery/results/freeze-d13/inputs/guccione_skin-ms-60mb-acgt-v2 guccione_skin-locate-500kb-v2 --limit 500000
