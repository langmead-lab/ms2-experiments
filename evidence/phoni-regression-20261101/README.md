# PHONI and Augmented PHONI regression check, Rockfish, 2026-11-01

This check shows that the PHONI and Augmented PHONI fixes are in the pinned freeze binaries
(`phoni` sha256 aeb7fb3b, `aug_phoni` sha256 79d50bf5). It ran in a lab scratch directory on
Rockfish; the files here are copied from it unchanged, so they name Rockfish paths.

- `setup.sh` cloned each tool and built the upstream base (PHONI 1b81e468, Augmented PHONI 74395806)
  and the pull-request head (2722db95, 0bf0b0c2) with spack GCC 13.1.0, `-include cstdint`, with and
  without `-ftrivial-auto-var-init=zero` (the pull-request heads only with it). `build.status` holds
  the exit status of each build.
- `setup.log` shows that the first checkout of the two pull-request heads failed (`reference is not
  a tree`), because the local clones did not yet have those commits. They were then checked out by
  hand at 2722db95 and 0bf0b0c2 and rebuilt; `run.log` records those heads.
- `run.sh` ran each pull request's `test/regression/check.py` (a 16-character reference and 7
  queries) with the query binary of each build, and with the freeze query binary in place of it on
  indexes built by the base build, the pull-request head build and the benchmark's own build.
- `run.log` holds the results. The freeze binaries answered 7 of 7 queries correctly on every index.
  The upstream base builds answered 4 of 7 (PHONI with zero initialization), 2 of 7 (PHONI without
  it) and 1 of 7 (Augmented PHONI, either way). The pull-request head builds answered 7 of 7.
- `objdump-compare.txt`, made on 2026-11-03 with `objcmp.sh`, compares each freeze binary with the
  pull-request head build of this check. The machine code is the same; the files differ only in the
  source paths that assertion macros embed, which are longer in this check's build directory.

On 2026-11-03, `ms2exp/build/build_legacy.sh --freeze-paths` rebuilt both tools from the upstream
base commits with `patches/phoni.patch` and `patches/aug_phoni.patch`; see `TOOLS.md` for the result.
