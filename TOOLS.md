# Tools

Every tool the paper's benchmarks ran, pinned in `tools.lock.json`. The freeze manifests
(`results/*/manifest.json`) pin each binary again by sha256 as an input, so a run cannot use any
other bytes. There are no git submodules here: `python -m ms2exp.build.tools --lock tools.lock.json`
reads the pins below in place of its flags.

## ms2

`langmead-lab/ms2` at commit 2885ee86, tagged `paper-v1`, with the submodules that commit records:
Orbit (`BenLangmead/runperm`) at bc4b455b and TeraTools (`BenLangmead/TeraTools`) at de922b03. There
are no local changes. `ms2exp.build.tools` builds three binaries from it with
`g++ -static-libstdc++ -static-libgcc`, so that one set of bytes runs on Rockfish and on Amazon Linux:

- `ms2-bench`: `make ms2-bench` in `src/`, assertions off. Every timed cell runs it.
- `ms2-stats`: the counting build, `ms2-bench`'s command plus `-DMS_STATS -DTMS_STATS`. The Fig ops
  operation counts come from it. The Makefile's own `ms2-stats` target sets only `-DMS_STATS`, so the
  build compiles the Makefile's sources itself.
- `ms2-replog`: `make ms2-replog`, the counting build plus `-DTMS_REP_LOG`, which logs one record per
  reposition for Fig rep.

It also compiles `ms2exp/build/ms_check.cpp` (`ms_check`), the whole-output MS checker. The freeze
built all four on AWS from git bundles of the three repositories, since the ms2 repository was
then private. The build ran on an i4i.metal box under Amazon Linux 2023, with its system compiler:
GCC 11.5.0 (Red Hat 11.5.0-5, packages `gcc-c++` and `libstdc++-static`), which each binary records
in its `.comment` section. Their sha256s are in `tools.lock.json` and in
`results/freeze-2885ee8-3c8c8967/build-tools/outputs.json`. The rerun reuses the same binaries.
The build's own record (`tools/build-2885ee8-ms2.json`) stayed on the box, since the build job
uploaded only the binaries; `ms2exp.build.job` uploads such records beside `outputs.json` for later
runs. On 2026-11-03 the same command, run in an Amazon Linux 2023 container at the freeze paths with
GCC 11.5.0-5 and the kernel headers of the box's kernel (6.1.189-236.400), produced all four binaries
byte for byte. Its record is `results/freeze-2885ee8-3c8c8967/build-tools/build-2885ee8-ms2.rebuilt.json`,
and `evidence/ms2-rebuild-20261103/` says how it ran. The kernel headers matter: with the 6.18
headers that the container's package resolver chose, the `PerfCounters` constructor of `ms2-bench`
differed, since `struct perf_event_attr` is larger there.

## ropebwt3

`lh3/ropebwt3` at b408339, built with plain `make -j2` (the Makefile's own flags: `-g -Wall
-Wc++-compat -O3` and OpenMP) using GCC 13.1.0, Rockfish's spack build, in
`/scratch16/blangme2/ms-bakery/ropebwt3` on 2026-09-21. Every `.c`, `.h` and `Makefile` of the
source the binary was built from matches upstream at that commit. The binary's sha256 is cfad8475.
It loads `libgomp.so.1`, which the manifests ship as `tools/runtime/libgomp.so.1` for hosts without it.

## MONI

`maxrossi91/moni` at 7f3f954. The changes are needed only to compile with GCC 13: `#include <array>`
at the top of `src/compress_dictionary.cpp` (an index construction helper, not part of the `ms`
query path), and `CMAKE_CXX_FLAGS=-include cstdint`. The pinned binary (`moni-ms`, sha256 4ea254dd)
is the `src/ms` of that build, built with spack GCC 13.1.0 and CMake 3.27.7 in Release mode
(`-O3 -DNDEBUG`). `patches/moni.patch` holds the include. MONI is GPL-3.0-or-later.

### Changes to fetched dependencies

MONI, PHONI and Augmented PHONI fetch about 15 libraries with CMake's FetchContent when they are
configured, several of them by branch name (`master`, `develop`, `main`), so a configure today can
fetch other code. `tools.lock.json` (`legacy_build.fetchcontent`) records the commit of every
dependency of the freeze builds. Two of those dependencies carried a one-line change, the same in
all three builds, needed only to compile with GCC 13:

- bwt2lcp (`nicolaprezza/bwt2lcp` at c08e6885): `#include <cstdint>` at the top of
  `internal/include.hpp` (`patches/bwt2lcp-cstdint.patch`);
- sux, a submodule of ShapedSlp (`koeppl/ShapedSlp` at 883f700e, sux at 73f24792): `#include <utility>`
  at the top of `sux/util/Vector.hpp` (`patches/sux-utility.patch`).

No other fetched dependency had a change to a tracked file.

## PHONI

`koeppl/phoni` at 1b81e46, with one fix: matching statistics start at the first occurrence of the
last character of the query. The query loops chose the second occurrence (a zero-based select called
with 1), so the BWT row and its sample could refer to different suffixes. The fix is upstream pull
request [koeppl/phoni#4](https://github.com/koeppl/phoni/pull/4) (head 2722db95), which adds a
regression test, and `patches/phoni.patch` holds the same commits. Built with
`CMAKE_CXX_FLAGS=-include cstdint -ftrivial-auto-var-init=zero` (spack GCC 13.1.0, Release). Without
the zero initialization the index builder reads uninitialized memory and can write a bad index,
which is a separate upstream bug. The pinned binary is `phoni`, sha256 aeb7fb3b. MIT license.

## Augmented PHONI

`drnatebrown/aug_phoni` at 7439580, with two fixes to the initialization of matching statistics:
the succeeding run and the initial occurrence. The fixes are upstream pull request
[drnatebrown/aug_phoni#2](https://github.com/drnatebrown/aug_phoni/pull/2) (head 0bf0b0c2), and
`patches/aug_phoni.patch` holds the same commits. The build flags are PHONI's. The pinned binary is
`aug_phoni`, sha256 79d50bf5. MIT license.

## Checks of the PHONI and AugPHONI binaries

The freeze binaries were built on September 22 from the local branches
`backup/pre-split-ms-initial-state` (9553550) and `backup/pre-split-aug-ms-initialization`
(d1b18d6). The pull request branches were split from those later, and their headers differ from
the branches the binaries were built from only in comments and whitespace.

To show that the fixes are in the pinned bytes, each pinned binary ran the regression test of its
pull request (`test/regression/check.py`: the 16-character reference ATTACCCTTGCTAAGT and 7
queries). The index was built with the upstream base build, the pull request build and the
benchmark's own build in turn. With each index the freeze binaries got 7 of 7 queries right. The
unpatched upstream builds got 4 of 7 (PHONI) and 1 of 7 (Augmented PHONI). The check ran on Rockfish
on 2026-11-01; its scripts and logs are in `evidence/phoni-regression-20261101/`.

The local branches are not public. The pull-request heads stand in for them: built on 2026-11-01
with the same compiler and flags, they give the same machine code as the freeze binaries. Their
`.text` sections have the same size and disassemble to the same instructions once addresses are
masked, and the files differ only in `.rodata`, which holds the source paths that assertion
macros embed and which is 64 bytes longer in the rebuild (`objdump-compare.txt` there).

## Rebuilding the legacy tools

`ms2exp/build/build_legacy.sh` builds ropebwt3, MONI, PHONI and Augmented PHONI from the commits,
patches and dependency pins in `tools.lock.json` (`--help` lists the options; its header shows the
Rockfish compiler paths). PHONI and Augmented PHONI are built from their upstream base commits with
the pull-request patches. `--freeze-paths` maps the build's directories onto the freeze build's
(`legacy_build.tools.*.freeze_paths`), so that the source paths the binaries embed come out the
same. Without it, the binaries differ from the pinned ones in those embedded paths, as
`evidence/phoni-regression-20261101/objdump-compare.txt` shows for PHONI and Augmented PHONI.

On 2026-11-03, `build_legacy.sh --freeze-paths` with spack GCC 13.1.0 and CMake 3.27.7 on
devlangmead1 produced all four pinned binaries byte for byte: ropebwt3 cfad8475, `moni-ms`
4ea254dd, `phoni` aeb7fb3b and `aug_phoni` 79d50bf5 (`evidence/legacy-rebuild-20261103/`). So the
bytes that were timed can be rebuilt from public commits and the patches here, and the pull-request
heads give the same binaries as the local branches the freeze builds came from.


## Preparation tools

The tools that prepared the inputs, but that the freeze manifests do not run, are pinned under
`prep_tools` in `tools.lock.json`, and `data/README.md` says which step used each. Three TeraLCP
builds made the index precursors, none of them the TeraTools commit that ms2 `paper-v1` pins
(de922b03): bce9ea7 the chr22-h10 and chr22-h100 minima (renumbered by `fmd2rlbwt`), 414f809 the
HPRC2 and mtb152 minima and phi and the chr22 phi, and af42783 (TeraLCP-minima) the Fig rep ladder
rungs, which were built from text with `data/prep/build_rlbwt.cpp` rather than from an FMD. Commits
bce9ea7 and 414f809 are not on GitHub yet (`prep_tools` records `public: false`).

## Runtime libraries

MONI, PHONI and Augmented PHONI link the spack GCC 13.1 runtime (`libstdc++.so.6`,
`libgcc_s.so.1`), which Amazon Linux 2023 lacks, and ropebwt3 loads that build's OpenMP runtime
(`libgomp.so.1`). The manifests ship Rockfish's copies as `tools/runtime/*`, byte-identical to the
files in the spack GCC 13.1.0 `lib64` directory, and the runner puts them on `LD_LIBRARY_PATH`. The
legacy binaries also carry that directory as their RPATH, so on Rockfish they need no
`LD_LIBRARY_PATH`. The GCC runtime libraries are GPL-3.0 with the GCC Runtime Library Exception.

All of these binaries, and the runtime libraries, are pinned by sha256 in `tools.lock.json` and in
the manifests (`ms2exp.build.tools --runtime` checks the runtime libraries against those pins); they are not yet public, so `python3 -m ms2exp.fetch --set toy` cannot obtain them outside the lab (`--dry-run` lists them). `build_legacy.sh
--static-runtime` builds MONI, PHONI and Augmented PHONI with the C++ runtime linked statically,
for hosts without a GCC 13 runtime; those binaries differ from the pinned ones by construction.
