#!/bin/bash
# Builds ropebwt3, MONI, PHONI and Augmented PHONI from the commits, patches and dependency pins in
# tools.lock.json. The options are those of ms2exp/build/legacy.py (run with --help).
#
# On Rockfish or devlangmead1, with the freeze builds' compiler and CMake:
#   G=/data/apps/extern/spack_on/gcc/9.5.0/gcc/13.1.0-nmocl3rvllh45gzlsritop3pidy63n2r/bin
#   C=/data/apps/extern/spack_on/gcc/9.3.0/cmake/3.27.7-iygsadyhjln3t7tr73pq56arnaxr3cgg/bin
#   PATH=$G:$PATH ms2exp/build/build_legacy.sh --cc $G/gcc --cxx $G/g++ --cmake $C/cmake \
#       --work DIR/work --out DIR/out --freeze-paths
set -euo pipefail
repo=$(cd "$(dirname "$0")/../.." && pwd)
PYTHONPATH="$repo${PYTHONPATH:+:$PYTHONPATH}" exec python3 -m ms2exp.build.legacy "$@"
