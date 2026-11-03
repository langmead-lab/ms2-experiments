#!/bin/bash
# The Fig ops query sets before their D13 sampling and caps (run_d13.sh prep_a): the held-out genomes as
# q-named workloads with truth, their 100 kbase pieces, and simulated PacBio HiFi reads from the same
# genomes with truth. Made on Rockfish Oct 25 to 29 2026; records in
# data/evidence/rockfish/results/campaign-20261012/inputs/{heldout,fig_ops,fig_ops2}/. Truth runs ropebwt3
# (sha256 cfad8475) on each reference's FMD and needs about 15 GB of memory for human472.
# Run from the repository root.
# Usage: data/prep/fig_ops_sets.sh
#   environment: INPUTS (destination, the campaign inputs/ directory), APR008 (GCA_050494745.1 genomic
#   .fna.gz), MTB152_PATTERNS (mtb152.5%.fa of data/prep/mtb152_split.sh), ROPEBWT3, HUMAN472_FMD,
#   MTB152_FMD (mtb152.95%.fmd), WORK, THREADS
set -euo pipefail
I=${INPUTS:?}; WORK=${WORK:?}; T=${THREADS:-16}
P() { python3 -m ms2exp.prep.pieces "$@"; }
S() { python3 -m ms2exp.sim.chr22_sim "$@"; }

# HPRC2: APR008 haplotype 1, 79 scaffolds split at N into 81 ACGT-only records.
H=$I/heldout/apr-APR008-hap1
mkdir -p $H
python3 -m ms2exp.prep.acgt_split "$APR008" $H/reads.fa $H/offsets.u64 $H/conversion_stats.json
P workload --reference human472 --name apr008-hap1-v1 --source $H/reads.fa --destination $I/fig_ops/human472-apr008-hap1-v1
P truth --workload $I/fig_ops/human472-apr008-hap1-v1 --ropebwt3 "$ROPEBWT3" --fmd "$HUMAN472_FMD" --work $WORK --threads $T
P pieces --workload $I/fig_ops/human472-apr008-hap1-v1 --length 100000 --destination $I/fig_ops/human472-apr008-hap1-p100k-v1

# mtb152: the 8 held-out genomes, each followed by its reverse complement.
P workload --reference mtb152 --name heldout16-v1 --source "$MTB152_PATTERNS" --destination $I/fig_ops/mtb152-heldout16-v1
P truth --workload $I/fig_ops/mtb152-heldout16-v1 --ropebwt3 "$ROPEBWT3" --fmd "$MTB152_FMD" --work $WORK --threads $T
P pieces --workload $I/fig_ops/mtb152-heldout16-v1 --length 100000 --destination $I/fig_ops/mtb152-heldout16-p100k-v1

# Simulated PacBio HiFi reads, the chr22 PacBio recipe of ms2exp.sim.chr22_sim on these sources.
S simulate --set pacbio --threads $T --cohort fig_ops2 --seed 2026102612 \
    --source $I/fig_ops/human472-apr008-hap1-v1/reads.fa --destination $I/fig_ops2/sim-human472-apr008-pacbio-60mb-v1 \
    --source-note "Fig Ops2 source: inputs/fig_ops/human472-apr008-hap1-v1/reads.fa, Arab Pangenome Reference APR008 haplotype 1 (GenBank GCA_050494745.1), absent from Human472, as 81 ACGT-only records"
S simulate --set pacbio --threads $T --cohort fig_ops2 --seed 2026102611 --pair-reverse-complements \
    --source $I/fig_ops/mtb152-heldout16-v1/reads.fa --destination $I/fig_ops2/sim-mtb152-pacbio-60mb-v1 \
    --source-note "Fig Ops2 source: inputs/fig_ops/mtb152-heldout16-v1/reads.fa, the 8 mtb152 genomes held out of the TeraTools 95 percent index (Nate Brown 95/5 split), each followed by its reverse complement; forward records only"
for r in human472-apr008 mtb152; do
    ref=${r%%-*}; fmd=$HUMAN472_FMD; [ $ref = mtb152 ] && fmd=$MTB152_FMD
    name=${r#*-}; from="APR008 haplotype 1"
    if [ $ref = mtb152 ]; then name=""; from="the forward records of the 8 held-out mtb152 genomes"; fi
    P workload --reference $ref --name ${name:+$name-}pacbio-60mb-v1 --cohort fig_ops2 \
        --source $I/fig_ops2/sim-$r-pacbio-60mb-v1/reads.fa --destination $I/fig_ops2/$r-pacbio-60mb-v1 \
        --source-note "Fig Ops2: simulated PacBio HiFi reads (chr22 PacBio recipe: pbsim3 3.0.5 multi-pass QSHMM-RSII, 10 passes, then ccs 6.4.0 defaults; template lengths from the real PacBio workload, about 60 Mbases of whole ACGT-only reads) from $from; read set and its manifest in inputs/fig_ops2/sim-$r-pacbio-60mb-v1"
    P truth --workload $I/fig_ops2/$r-pacbio-60mb-v1 --ropebwt3 "$ROPEBWT3" --fmd "$fmd" --work $WORK --threads $T
done
