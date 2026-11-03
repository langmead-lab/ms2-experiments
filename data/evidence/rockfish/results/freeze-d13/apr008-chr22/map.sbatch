#!/bin/bash
#SBATCH -J orbitms-apr008-chr22map
#SBATCH -p shared
#SBATCH -c 16
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH -o %x-%j.log
# Align every APR008 hap1 scaffold to CHM13 v2.0 chr22 to find the scaffold(s) that are its chr22.
set -euo pipefail
D=/scratch16/blangme2/ms-bakery/results/freeze-d13/apr008-chr22
M=/scratch16/blangme2/ms-bakery/tools/minimap2/minimap2
S=/scratch16/blangme2/ms-bakery/tools/hifi-env/.pixi/envs/default/bin/samtools
cd $D
$S faidx /data/blangme2/fasta/chm13_v2.0/chm13v2.0.fasta chr22 > chm13-chr22.fa
/usr/bin/time -v $M -x asm5 -t 16 --secondary=no chm13-chr22.fa \
  /scratch16/blangme2/ms-bakery/results/campaign-20261012/inputs/heldout/apr-APR008-hap1/raw/GCA_050494745.1_apr008.1_v1_genomic.fna.gz \
  > apr008-hap1-vs-chm13-chr22.paf 2> minimap2.stderr
# Per query scaffold: aligned bases on chr22 (mapq >= 20), and the target span it covers.
awk '$12>=20 {q[$1]=$2; a[$1]+=$4-$3; if(!($1 in lo)||$8<lo[$1])lo[$1]=$8; if($9>hi[$1])hi[$1]=$9}
     END {for(k in a) printf "%s\t%d\t%d\t%d\t%d\n", k, q[k], a[k], lo[k], hi[k]}' apr008-hap1-vs-chm13-chr22.paf \
  | sort -k3,3nr > per-scaffold.tsv
