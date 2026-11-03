# Patches

The source changes the freeze binaries of the legacy tools were built with. `ms2exp/build/build_legacy.sh`
applies all of them; `tools.lock.json` says where each one goes.

- `phoni.patch` and `aug_phoni.patch`: the PHONI and Augmented PHONI changes, as `git format-patch`
  output against their upstream base commits (`revisions.json`). They are the commits of the upstream
  pull requests koeppl/phoni#4 and drnatebrown/aug_phoni#2, kept here in case those are never merged.
  Both projects are MIT (copyright 2020 Massimiliano Rossi). Apply one with `git am` on a checkout of
  its base commit.
- `moni.patch`: `#include <array>` at the top of MONI's `src/compress_dictionary.cpp` (an index
  construction helper), needed to compile with GCC 13. MONI is GPL-3.0-or-later. Apply with
  `git apply` on a checkout of 7f3f954.
- `bwt2lcp-cstdint.patch` and `sux-utility.patch`: one-line includes, needed to compile with GCC 13,
  in two dependencies that CMake fetches when MONI, PHONI and Augmented PHONI are configured: bwt2lcp
  (`internal/include.hpp`, MIT, Nicola Prezza) and sux, a submodule of ShapedSlp
  (`vigna/sux`, `sux/util/Vector.hpp`, copyright Stefano Marchini, LGPL-3.0-or-later or GPL-3.0-or-later
  with the GCC Runtime Library Exception). Apply with `git apply` in the
  dependency's checkout (sux: in `external/sux` of ShapedSlp).
