# Rebuild of the legacy tools, devlangmead1, 2026-11-03

`run_legacy.sh` ran `ms2exp/build/build_legacy.sh --freeze-paths --jobs 1` with the freeze builds'
compiler and CMake (spack GCC 13.1.0, CMake 3.27.7) in a new scratch directory on devlangmead1,
from the repository's patches and `tools.lock.json`. It cloned ropebwt3, MONI, PHONI and Augmented
PHONI at their upstream commits from GitHub, applied `patches/moni.patch`, `patches/phoni.patch` and
`patches/aug_phoni.patch` (PHONI and Augmented PHONI therefore at their pull-request heads, not the
local branches the freeze binaries came from), cloned every fetched dependency at its pinned commit
(from GitHub and GitLab) with `patches/bwt2lcp-cstdint.patch` and `patches/sux-utility.patch`, and built.

All four binaries are byte-identical to the pinned freeze binaries (`result.txt`; `build-legacy.json`
is the script's record, with every command). It took 15 minutes with one build job, at a peak of
0.6 GB.
