#!/bin/bash
# Find the APR008 haplotype 1 scaffold that is chromosome 22: align every scaffold of GCA_050494745.1 to
# CHM13 v2.0 chr22 with minimap2 2.28 (-x asm5 --secondary=no) and summarize the alignments per scaffold.
# per-scaffold.tsv sums every alignment of mapping quality 20 or more (scaffold, length, aligned query
# bases, first and last chr22 coordinate covered). long-arm.tsv keeps only alignment blocks of at least
# 500 kbases; on the frozen run it named JBHMUJ010000002.1 alone over chr22 17.47 to 51.32 Mb, with
# 33,878,807 aligned bases. The frozen outputs are in data/evidence/rockfish/results/freeze-d13/apr008-chr22/.
# Usage: data/prep/apr008_chr22_map.sh OUTDIR CHM13_FASTA APR008_FNA_GZ   (environment: MINIMAP2, SAMTOOLS, THREADS)
set -euo pipefail
D=$1; CHM13=$2; APR008=$3
M=${MINIMAP2:-minimap2}
S=${SAMTOOLS:-samtools}
mkdir -p "$D"
cd "$D"
$S faidx "$CHM13" chr22 > chm13-chr22.fa
$M -x asm5 -t ${THREADS:-16} --secondary=no chm13-chr22.fa "$APR008" > apr008-hap1-vs-chm13-chr22.paf 2> minimap2.stderr
awk '$12>=20 {q[$1]=$2; a[$1]+=$4-$3; if(!($1 in lo)||$8<lo[$1])lo[$1]=$8; if($9>hi[$1])hi[$1]=$9}
     END {for(k in a) printf "%s\t%d\t%d\t%d\t%d\n", k, q[k], a[k], lo[k], hi[k]}' apr008-hap1-vs-chm13-chr22.paf \
  | sort -k3,3nr > per-scaffold.tsv
{
  printf 'scaffold\tscaffold_length\talignments\taligned_query_bases\tchm13_chr22_start\tchm13_chr22_end\n'
  awk '$12>=20 && $11>=500000 {q[$1]=$2; a[$1]+=$4-$3; n[$1]++; if(!($1 in lo)||$8<lo[$1])lo[$1]=$8; if($9>hi[$1])hi[$1]=$9}
       END {for(k in a) printf "%s\t%d\t%d\t%d\t%d\t%d\n", k, q[k], n[k], a[k], lo[k], hi[k]}' apr008-hap1-vs-chm13-chr22.paf \
    | sort -k4,4nr
} > long-arm.tsv
head -3 long-arm.tsv
