#!/bin/bash
# Rebuild of the legacy tools with the freeze compiler and path maps, under a 1.4 GB RSS watchdog.
D=/scratch16/blangme2/ms-bakery/scratch/legacy-rebuild-20261103
G=/data/apps/extern/spack_on/gcc/9.5.0/gcc/13.1.0-nmocl3rvllh45gzlsritop3pidy63n2r/bin
C=/data/apps/extern/spack_on/gcc/9.3.0/cmake/3.27.7-iygsadyhjln3t7tr73pq56arnaxr3cgg/bin
cd $D
export PATH=$G:$C:$PATH
setsid /usr/bin/time -v $D/repo/ms2exp/build/build_legacy.sh --cc $G/gcc --cxx $G/g++ --cmake $C/cmake \
  --work $D/work --out $D/out --freeze-paths --jobs 1 > $D/build.log 2>&1 &
pid=$!
pgid=$(ps -o pgid= $pid | tr -d ' ')
peak=0
while kill -0 $pid 2>/dev/null; do
  rss=$(ps -g $pgid -o rss= 2>/dev/null | awk '{s+=$1} END {print s+0}')
  [ "$rss" -gt "$peak" ] && peak=$rss
  if [ "$rss" -gt 1400000 ]; then echo "watchdog: rss $rss KiB, stopping" >> $D/watchdog.log; kill -- -$pgid; fi
  echo "$peak" > $D/peak_rss_kib
  sleep 5
done
wait $pid; echo "exit $?" >> $D/build.status
