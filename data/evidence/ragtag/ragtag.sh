#!/bin/bash -l
#SBATCH --job-name="{0}"
#SBATCH --output="{0}.out"
#SBATCH -p parallel
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=48
#SBATCH -A blangme2
#SBATCH --time 1:00:00

conda activate ragtag
cd /home/USER/data-blangme2/fasta/hprc-v2/final_update/ragtag/{0}
cp /data/blangme2/fasta/hprc-v2/final_update/fastas/{1} ./{0}.fa
/usr/bin/time -v ragtag.py scaffold /home/USER/data-blangme2/fasta/hprc-v2/final_update/fastas/100000_CHM13.pri.fa {0}.fa
