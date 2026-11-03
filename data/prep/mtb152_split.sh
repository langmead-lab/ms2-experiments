#!/bin/bash
# Reconstructed, not the original script. The mtb152 95/5 split and its RLBWT precursors.
#
# Nate Brown made the split on Oct 1 2026 with a setup script that is not in this repository. Its log,
# data/evidence/mtb152/setup.20261001_132202.log, is the only record, and this script repeats the commands
# that log shows. It has not been run; the sha256 checks at the end say whether it reproduces the frozen
# files. From Li's ropebwt3 index of the 152 M. tuberculosis genomes (mtb152.fmd, whose sequences are each
# genome and then its reverse complement):
#   patterns: ropebwt3 get sequences 0 to 15 (the first 8 genomes and their reverse complements), then N
#     characters removed (the log says "Cleaning pattern of N characters"; how is not recorded, and this
#     removes them from the sequence lines), giving mtb152.5%.fa (16 records with ropebwt3 get headers such
#     as ">0 1318782311", one sequence line each);
#   index: ropebwt3 get sequences 16 to 302 in steps of 2 (the forward strands of the other 144 genomes),
#     then ropebwt3 build -t 32 -do mtb152.95%.fmd (both strands are added back);
#   precursors: fmd2rlbwt and TeraLCP -f fmd -ominima -ophi on that FMD, the dump_runs and teralcp rules of the
#     terams-fig Snakemake workflow (TeraLCP 414f809 there).
# The log does not record where mtb152.tar.gz came from beyond Li's ropebwt3 index collection; MTB152_TAR names
# the downloaded archive.
# Usage: data/prep/mtb152_split.sh OUTDIR   (environment: MTB152_TAR, ROPEBWT3, FMD2RLBWT, TERALCP, THREADS)
set -euo pipefail
OUT=$1; T=${THREADS:-32}; R=${ROPEBWT3:-ropebwt3}
mkdir -p $OUT/pattern $OUT/reference $OUT/orbit $OUT/tmp
cd $OUT
tar xzf "$MTB152_TAR"
# Nate's copy of mtb152.fmd (dated Jun 15 2024 beside his log) has this sha256.
test -e mtb152.fmd || $R build -i mtb152.fmr -do mtb152.fmd
check() {  # file expected-sha256
    if [ "$(sha256sum "$1" | cut -d' ' -f1)" = "$2" ]; then echo "$1 equals the frozen file"; else echo "$1 differs from the frozen file"; fi
}
check mtb152.fmd 974f61b02e3e5aed5e6d394834440a5cf9ad9cf920026b7c93c85829110c81ed
test "$($R stat mtb152.fmd | awk '/sequences/ {print $1}')" = 304
for i in $(seq 0 15); do $R get mtb152.fmd $i > pattern/$i.fa; done
for i in $(seq 0 15); do cat pattern/$i.fa; done | awk '/^>/ {print; next} {gsub(/[Nn]/, ""); print}' > pattern/mtb152.5%.fa
for i in $(seq 16 2 302); do $R get mtb152.fmd $i; done > reference/mtb152.95%.fa
$R build -t $T -do mtb152.95%.fmd reference/mtb152.95%.fa
rm -f reference/mtb152.95%.fa
$FMD2RLBWT mtb152.95%.fmd orbit/mtb152 > orbit/mtb152.fmd2rlbwt.txt
$TERALCP -f fmd -i mtb152.95%.fmd -t tmp/t -p 16 -ominima orbit/mtb152.min -ophi orbit/mtb152.phi
check pattern/mtb152.5%.fa 7df00bdd032400d635d99bb1a933dcd8d03eff1d51471def15b41fa6399dc659
check mtb152.95%.fmd 181c92baeda3f414e50c1e5e87b71ac2c3e54d8bc47237da151aecc192b03381
check orbit/mtb152.bwt.heads 7eaaa9ae908a1bde8cd6c79de8a862f2828513a0220af951420682a4d6c09216
check orbit/mtb152.bwt.len 024e9f54f656eb08fe2b36927006edcc25023fc677653ba78d34641af8f13d88
check orbit/mtb152.min 61431b66c9ebc1e6a8e8006917aa7062b08fd7d93ee59683037c6f7cea84caf7
check orbit/mtb152.phi 6598a9b5a38438f84693cefc1f5a37d6136cb5e921903e354f98263f7ee8fa7c
