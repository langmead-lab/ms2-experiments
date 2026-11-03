#!/bin/bash
# Fresh clones and builds for the freeze-binary regression check.
set -u
W=/scratch16/blangme2/ms-bakery/scratch/ms2exp-binary-check
U=/scratch16/blangme2/ms-bakery/scratch/upstream-pr
mkdir -p $W/logs $W/build && cd $W
for r in phoni aug_phoni; do
  [ -d $r ] || git clone -q $U/$r $r
  git -C $r fetch -q $U/$r 'refs/heads/*:refs/remotes/src/*'
done
git -C phoni rev-parse src/pr/ms-initial-state src/phoni
git -C aug_phoni rev-parse src/pr/aug-ms-initialization src/phoni
# Separate worktree-free source copies per revision
mk() { # repo rev name
  [ -d src-$3 ] || { git clone -q $1 src-$3 && git -C src-$3 checkout -q $2; }
  echo "src-$3 at $(git -C src-$3 rev-parse HEAD)"
}
mk phoni 1b81e468086f99a1428f8af190835659bf47f503 phoni-base
mk phoni 2722db95d379412075495996d51cd0671a7e97e7 phoni-pr
mk aug_phoni 74395806673e3ddc4d83ac5261133d64c4885ba3 aug_phoni-base
mk aug_phoni 0bf0b0c2cae1aa68560e25f2973fc0b91deb1d59 aug_phoni-pr
common="newscanNT.x pscan.x newscan.x pfp_thresholds pfp_thresholds64 compress_dictionary rlbwt bwtparse pfbwtNT.x"
for s in phoni-base phoni-pr aug_phoni-base aug_phoni-pr; do
  case $s in phoni*) t="$common phoni build_phoni";; *) t="$common aug_phoni build_aug pfp_lcp";; esac
  for v in zero plain; do
    [ $s = ${s%-pr} ] || [ $v = zero ] || continue   # PR head only with zero init
    f="-include cstdint"; [ $v = zero ] && f="$f -ftrivial-auto-var-init=zero"
    (cmake -S src-$s -B build/$s-$v -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_FLAGS="$f" > logs/$s-$v.cfg 2>&1 && \
     cmake --build build/$s-$v -j4 --target $t > logs/$s-$v.build 2>&1; echo "$s $v exit=$?" >> logs/build.status) &
  done
done
wait
cat logs/build.status
