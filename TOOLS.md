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
built all four on AWS from git bundles of the three repositories, since the ms2 repository is
private. Their sha256s are in `tools.lock.json` and in
`results/freeze-2885ee8-3c8c8967/build-tools/outputs.json`. The rerun reuses the same binaries.

## ropebwt3

`lh3/ropebwt3` at b408339, built with plain `make -j2`. Every `.c`, `.h` and `Makefile` of the
source the binary was built from matches upstream at that commit. The binary's sha256 is cfad8475.
It loads `libgomp.so.1`, which the manifests ship as `tools/runtime/libgomp.so.1` for hosts without it.

## MONI

`maxrossi91/moni` at 7f3f954. The changes are needed only to compile with GCC 13: `#include <array>`
at the top of `src/compress_dictionary.cpp` (an index construction helper, not part of the `ms`
query path), and `CMAKE_CXX_FLAGS=-include cstdint`. The pinned binary (`moni-ms`, sha256 4ea254dd)
is the `src/ms` of that build, built with spack GCC 13.1 in Release mode. MONI is GPL-3.0-or-later.

## PHONI

`koeppl/phoni` at 1b81e46, with one fix: matching statistics start at the first occurrence of the
last character of the query. The query loops chose the second occurrence (a zero-based select called
with 1), so the BWT row and its sample could refer to different suffixes. The fix is upstream pull
request [koeppl/phoni#4](https://github.com/koeppl/phoni/pull/4) (head 2722db95), which adds a
regression test, and `patches/phoni.patch` holds the same commits. Built with
`CMAKE_CXX_FLAGS=-include cstdint -ftrivial-auto-var-init=zero` (spack GCC 13.1, Release). Without
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
on 2026-11-01, in `/scratch16/blangme2/ms-bakery/scratch/ms2exp-binary-check` (lab-internal).

## Runtime libraries

MONI, PHONI and Augmented PHONI link the spack GCC 13.1 runtime (`libstdc++.so.6`,
`libgcc_s.so.1`), which Amazon Linux 2023 lacks. The manifests ship Rockfish's copies as
`tools/runtime/*` and the runner puts them on `LD_LIBRARY_PATH`. Their sha256s are in
`tools.lock.json`.
