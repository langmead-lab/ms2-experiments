#!/bin/bash
# Run each PR's committed regression test (test/regression/check.py at the PR head) against
# fresh base builds, fresh PR-head builds, and the pinned freeze query binaries.
W=/scratch16/blangme2/ms-bakery/scratch/ms2exp-binary-check
P=/scratch16/blangme2/ms-bakery/.venv/bin/python
FZ=/scratch16/blangme2/ms-bakery/results/real-hprc2-1d580b6/chr22/bin
cd $W
# mirror <index-pipeline build dir> <name> <query exe name> <query binary>
mirror() {
  d=build/$2; rm -rf $d; mkdir -p $d/test/src
  for f in $1/*; do b=$(basename $f); [ $b = test ] || ln -s $(readlink -f $f) $d/$b; done
  for f in $1/test/src/*; do b=$(basename $f); [ $b = $3 ] || ln -s $(readlink -f $f) $d/test/src/$b; done
  ln -s $4 $d/test/src/$3
}
mirror build/phoni-base-zero phoni-freeze-on-base-zero phoni $FZ/phoni
mirror build/phoni-pr-zero phoni-freeze-on-pr-zero phoni $FZ/phoni
mirror /scratch16/blangme2/ms-bakery/phoni/build-bakery phoni-freeze-on-bakery phoni $FZ/phoni
mirror build/aug_phoni-base-zero aug_phoni-freeze-on-base-zero aug_phoni $FZ/aug_phoni
mirror build/aug_phoni-pr-zero aug_phoni-freeze-on-pr-zero aug_phoni $FZ/aug_phoni
mirror /scratch16/blangme2/ms-bakery/aug_phoni/build-bakery aug_phoni-freeze-on-bakery aug_phoni $FZ/aug_phoni
sha256sum $FZ/phoni $FZ/aug_phoni
for r in phoni aug_phoni; do
  T=src-$r-pr/test/regression/check.py
  echo "test: $T at $(git -C src-$r-pr rev-parse HEAD), check.py sha256 $(sha256sum < $T | cut -c1-16)"
  for b in $r-base-zero $r-base-plain $r-pr-zero $r-freeze-on-base-zero $r-freeze-on-pr-zero $r-freeze-on-bakery; do
    q=build/$b/test/src/$r
    echo "== $r $b query=$(readlink -f $q) sha256=$(sha256sum < $(readlink -f $q) | cut -c1-16)"
    $P $T build/$b 2>/dev/null; echo "exit=$?"
  done
done
