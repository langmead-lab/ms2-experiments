#!/bin/bash
# The A1 mtb152 error series (mtb152-A2-*): 20 Mbase batches from the 8 forward held-out mtb152 genomes,
# error-free or with independent per-base errors injected by ms2exp.figures.fig_rep simulate, with the
# seeds of the frozen sets (manifests in data/evidence/rockfish/results/campaign-20261012/fig_rep/inputs/).
# The rates make the mean psi steps per reposition fall by about the same factor from one batch to the next.
# Only the read sets are made here; the freeze manifest runs them. Run from the repository root.
# Usage: data/prep/series_a2.sh SOURCE_READS_FA DESTINATION_DIR
#   SOURCE_READS_FA is fig_ops/mtb152-heldout16-v1/reads.fa of data/prep/fig_ops_sets.sh.
set -euo pipefail
SRC=$1; W=$2
batch() {  # name method rate seed
    local name=mtb152-A2-$1
    test -e $W/$name/manifest.json && return 0
    python3 -m ms2exp.figures.fig_rep simulate --source $SRC --pair-reverse-complements \
        --source-note "8 forward held-out mtb152 genomes" --destination $W/$name --method $2 --rate $3 \
        --bases 20000000 --seed $4
}
batch exact exact 0 2026102801
batch e0.000032 inject 0.000032 2026102802
batch e0.00011 inject 0.00011 2026102803
batch e0.00019 inject 0.00019 2026102804
batch e0.00030 inject 0.00030 2026102805
batch e0.00047 inject 0.00047 2026102806
batch e0.00076 inject 0.00076 2026102807
batch e0.0014 inject 0.0014 2026102808
batch e0.0030 inject 0.0030 2026102809
batch e0.0050 inject 0.0050 2026102810
