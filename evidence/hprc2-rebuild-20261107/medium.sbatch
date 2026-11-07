#!/bin/bash
#SBATCH -J orbitms-hprc2-rebuild-medium
#SBATCH -p shared,parallel
#SBATCH -A blangme2
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH -o %x-%j.out
# Medium tier: TeraLCP built from the public tags reproduces the deposited mtb152 and chr22 precursors.
#   mtb152: TeraLCP ms2-precursors (414f809) -f fmd -ominima -ophi on the deposited mtb152 FMD, as in
#     data/prep/mtb152_split.sh; fmd2rlbwt for the run files.
#   chr22-h10, chr22-h100: data/prep/chr22_precursors.sh with TeraLCP ms2-chr22-minima (bce9ea7) for the
#     minima and ms2-precursors for the phi samples, on the deposited chr22 FMDs.
# Usage: sbatch medium.sbatch REPO [TERATOOLS_REPO]
set -euo pipefail
REPO=$1; export TERATOOLS_REPO=${2:-https://github.com/BenLangmead/TeraTools.git}
W=/tmp/orbitms-hprc2-medium-$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
df -h /tmp; hostname; gcc --version | head -1
export MS2EXP_DATA=$W/data THREADS=16
D=$MS2EXP_DATA
cd $REPO
python3 -m ms2exp.fetch --data $D tools/ropebwt3-cfad8475/ropebwt3 tools/runtime/libgcc_s.so.1 tools/runtime/libgomp.so.1 \
    ropebwt3/mtb152/index.fmd ropebwt3/chr22-h10/index.fmd ropebwt3/chr22-h100/index.fmd
chmod +x $D/tools/ropebwt3-cfad8475/ropebwt3
data/prep/rebuild_hprc2.sh $W/work tools
TL414=$W/work/teratools/src/TeraLCP/TeraLCP; F2R=$W/work/runperm/src/ms/tools/fmd2rlbwt/fmd2rlbwt
# The same TeraTools tree checked out at the chr22 minima tag, with its own build.
git clone -q $W/work/teratools $W/work/teratools-bce9ea7
git -C $W/work/teratools-bce9ea7 checkout -q ms2-chr22-minima
test "$(git -C $W/work/teratools-bce9ea7 rev-parse HEAD)" = bce9ea7f387b71dd413beb15b8236cd705aa83b0
cp -r $W/work/teratools/ropebwt3 $W/work/teratools-bce9ea7/ropebwt3
make -C $W/work/teratools-bce9ea7/src thirdparty
make -C $W/work/teratools-bce9ea7/src/TeraLCP TeraLCP
TLBCE=$W/work/teratools-bce9ea7/src/TeraLCP/TeraLCP
sha256sum $TL414 $TLBCE $F2R | tee $SLURM_SUBMIT_DIR/medium-binaries.sha256
want() { awk -F'\t' -v n="$1" '$1 == n {print $3}' data/inputs.tsv; }
report() {  # NAME FILE
    local got; got=$(sha256sum "$2" | cut -d' ' -f1)
    if [ "$got" = "$(want $1)" ]; then echo "equal   $1 $got"; else echo "DIFFER  $1 $got expected $(want $1)"; fi
}
R=$SLURM_SUBMIT_DIR/medium-result.txt; : > $R
M=$W/mtb152; mkdir -p $M/tmp
$F2R $D/ropebwt3/mtb152/index.fmd $M/mtb152 > $M/fmd2rlbwt.txt
/usr/bin/time -v $TL414 -f fmd -i $D/ropebwt3/mtb152/index.fmd -t $M/tmp/t -p 16 -ominima $M/mtb152.min -ophi $M/mtb152.phi \
    > $M/teralcp.stdout 2> $SLURM_SUBMIT_DIR/medium-mtb152-teralcp.time
for f in bwt.heads bwt.len; do report precursors/mtb152/$f $M/mtb152.$f | tee -a $R; done
report precursors/mtb152/min $M/mtb152.min | tee -a $R
report precursors/mtb152/phi $M/mtb152.phi | tee -a $R
for c in chr22-h10 chr22-h100; do
  FMD2RLBWT=$F2R TERALCP_MINIMA=$TLBCE TERALCP_PHI=$TL414 ROPEBWT3=$D/tools/ropebwt3-cfad8475/ropebwt3 \
    LD_LIBRARY_PATH=$D/tools/runtime /usr/bin/time -v data/prep/chr22_precursors.sh $D/ropebwt3/$c/index.fmd $W/$c \
    2> $SLURM_SUBMIT_DIR/medium-$c.time
  for f in bwt.heads bwt.len min phi; do report precursors/$c/$f $W/$c/fmd.$f | tee -a $R; done
done
