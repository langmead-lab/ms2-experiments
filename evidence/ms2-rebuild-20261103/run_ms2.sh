#!/bin/bash
# ms2 freeze build again with kernel-headers 6.1.189-236.400, in an Amazon Linux 2023 container at the freeze paths, under a 1.4 GB RSS watchdog.
D=/scratch16/blangme2/ms-bakery/scratch/legacy-rebuild-20261103/ms2
cd $D
I=/work/stage/inputs
setsid singularity exec --cleanenv --no-home --bind $D/bind2/work:/work --bind $D/tmp:/tmp $D/al2023k \
  bash -c "g++ --version | head -1; python3 /work/code/bench/freeze/build_tools.py --parts ms2 --ms2-repo $I/cc88d25eec54d2cdb21f2e934cafba53fc4b4c55a24af722dc098703d1bf323e/ms2.bundle --submodule-repo orbit $I/f13bff97b6f373b18da426da9130966a8c9ad380eee91a3b3d4dfadae592a960/orbit.bundle --submodule-repo teratools $I/c941f5cf0ce8108d0087c9a826bad138c45b688d82d30e70356dfd97ee986ed1/teratools.bundle --ms2-commit 2885ee86270b4187fc45dcc69efd7de48a4aac03 --root /work/stage/work/freeze-2885ee8-3c8c8967/build-tools/outputs --work /work/stage/work/freeze-2885ee8-3c8c8967/build-tools/scratch/ms2-src" > $D/build2.log 2>&1 &
pid=$!
pgid=$(ps -o pgid= $pid | tr -d ' ')
peak=0
while kill -0 $pid 2>/dev/null; do
  rss=$(ps -g $pgid -o rss= 2>/dev/null | awk '{s+=$1} END {print s+0}')
  [ "$rss" -gt "$peak" ] && peak=$rss
  if [ "$rss" -gt 1400000 ]; then echo "watchdog: rss $rss KiB, stopping" >> $D/watchdog2.log; kill -- -$pgid; fi
  echo "$peak" > $D/peak_rss_kib2
  sleep 2
done
wait $pid; echo "exit $?" >> $D/build2.status
