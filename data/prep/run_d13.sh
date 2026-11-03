#!/bin/bash
# The D13 read sets of the freeze (BENCH_FREEZE_DESIGN.md D13), in the order they were made on Rockfish
# on Oct 30 2026 by the three Slurm jobs in data/evidence/rockfish/results/freeze-d13/prep-*.sbatch.
# Each step is a ms2exp.build.d13_prep subcommand; a step whose destination exists stops, so run this on
# an empty --out. Sources are the paths d13_prep.py names: the pools of ms2exp.prep.pool, the
# explore-acgt-60mb-v1 sets of ms2exp.prep.acgt_workload, the Fig ops sets of data/prep/fig_ops_sets.sh,
# the APR008 assembly and the ERR14066769 download. Run from the repository root.
# Usage: data/prep/run_d13.sh [all|a|chr22|skin]   (environment: OUT, CAMPAIGN, THREADS)
set -euo pipefail
OUT=${OUT:-/scratch16/blangme2/ms-bakery/results/freeze-d13/inputs}
CAMPAIGN=${CAMPAIGN:-/scratch16/blangme2/ms-bakery/results/campaign-20261012}
THREADS=${THREADS:-16}
X() { python3 -m ms2exp.build.d13_prep --out "$OUT" "$@"; }
I=$CAMPAIGN/inputs

prep_a() {
    # Sets that slice existing truth: 60 Mb caps, the Fig ops genome samples, the real short set, then
    # the locate prefixes.
    X cap $I/human_pacbio/explore-acgt-60mb-v1 human_pacbio-ms-60mb-acgt-v2
    X cap $I/human_ont_modern/explore-acgt-60mb-v1 human_ont_modern-ms-60mb-acgt-v2
    X cap $I/fig_ops2/human472-apr008-pacbio-60mb-v1 human472-apr008-pacbio-60mb-v2
    X cap $I/fig_ops2/mtb152-pacbio-60mb-v1 mtb152-pacbio-60mb-v2
    X sample $I/fig_ops/human472-apr008-hap1-p100k-v1 human472-apr008-hap1-p100k-s600-v1 --seed 2026103001
    X sample $I/fig_ops/mtb152-heldout16-p100k-v1 mtb152-heldout16-p100k-s600-v1 --seed 2026103002
    X short-60mb
    X cap $OUT/human_short-ms-60mb-acgt-v2 human_short-locate-6mb-v2 --limit 6000000
    X cap $OUT/human_pacbio-ms-60mb-acgt-v2 human_pacbio-locate-20mb-v2 --limit 20000000
    X cap $OUT/human_ont_modern-ms-60mb-acgt-v2 human_ont_modern-locate-10mb-v2 --limit 10000000
}

prep_chr22() {
    # The APR008 hap1 chr22 long arm (data/prep/apr008_chr22_map.sh identifies the scaffold), then simulated
    # short, PacBio and ONT reads with truth on chr22-h10 and chr22-h100.
    test -e $OUT/haplotype-APR008-hap1-chr22-v1 || X apr008-chr22 --scaffold JBHMUJ010000002.1 --evidence \
        "minimap2 2.28 -x asm5 of every APR008 hap1 scaffold against CHM13 v2.0 chr22 (job 32742987, freeze-d13/apr008-chr22): JBHMUJ010000002.1 covers chr22 17.4-51.3 Mb (33.85 Mb in alignments of at least 500 kb, 99.5 percent of the scaffold) and no other scaffold aligns there; chr22 0-17.4 Mb (short arm and centromere) is spread over at least 10 scaffolds shared with other acrocentrics and is not included"
    for s in short pacbio ont; do X simulate --set $s --threads $THREADS; done
}

prep_skin() {
    # Filter first, cap at 60,000,000 bases, new truth on the human472 FMD; then the 0.5 Mb locate prefix.
    X skin-60mb --threads $THREADS
    X cap $OUT/guccione_skin-ms-60mb-acgt-v2 guccione_skin-locate-500kb-v2 --limit 500000
}

case "${1:-all}" in
    a) prep_a ;;
    chr22) prep_chr22 ;;
    skin) prep_skin ;;
    all) prep_a; prep_chr22; prep_skin ;;
    *) echo "usage: $0 [all|a|chr22|skin]" >&2; exit 2 ;;
esac
