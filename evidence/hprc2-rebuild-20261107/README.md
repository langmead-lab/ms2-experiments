# Rebuild of the HPRC2 precursor steps from public code, Rockfish, 2026-11-07

`data/prep/rebuild_hprc2.sh` rebuilds the HPRC2 inputs that are not deposited. These runs checked, in two
tiers, that its steps reproduce the frozen files byte for byte. Both jobs ran on Rockfish compute nodes with
Spack GCC 13.1.0 and wrote everything large to node-local `/tmp`.

TeraTools came from a git bundle of the tags `ms2-precursors` (414f809) and `ms2-chr22-minima` (bce9ea7),
made before the tags were pushed to `BenLangmead/TeraTools` (`TERATOOLS_REPO` pointed at the bundle). A
bundle holds the same commits, so the checkout is the same tree a clone from GitHub gives. Everything else
came from its public source: fmd2rlbwt from `BenLangmead/runperm` 2a3bf10, ropebwt3 b408339 from
`lh3/ropebwt3` (in TeraTools' `ropebwt3/` directory), sdsl-lite from `simongog/sdsl-lite` (cloned at its
head by TeraTools' Makefile), the pinned ropebwt3 binary and runtime libraries from the `tools-v1` release
via `ms2exp.fetch`, the FMDs from deposit 2 (Zenodo 23964220) and `human472.fmr.gz` from Li's record
(Zenodo 14854401).

## Medium tier: deposited mtb152 and chr22 precursors (`medium.sbatch`, job 32986702)

TeraLCP built from `ms2-precursors` ran `-f fmd -ominima -ophi` on the deposited mtb152 FMD, as in
`data/prep/mtb152_split.sh`. `data/prep/chr22_precursors.sh` ran on the deposited chr22-h10 and chr22-h100
FMDs with TeraLCP built from `ms2-chr22-minima` for the minima and from `ms2-precursors` for the phi
samples. fmd2rlbwt made the run files in both cases.

All 12 files equal the deposited ones (`medium-result.txt`): `bwt.heads`, `bwt.len`, `min` and `phi` of
mtb152, chr22-h10 and chr22-h100.

| Step | Wall | Peak RSS |
|---|---|---|
| TeraLCP on mtb152, 16 threads | 54 s | 0.2 GB |
| `chr22_precursors.sh` chr22-h10, 16 threads | 3 min 51 s | 1.6 GB |
| `chr22_precursors.sh` chr22-h100, 16 threads | 18 min 27 s | 2.1 GB |

The rebuilt binaries (`medium-binaries.sha256`) differ from the sha256s in `tools.lock.json`, which belong
to the binaries built on 2026-10-01 in another build directory. Their outputs are identical.

## Cheap tier: HPRC2 FMD and run files (`cheap.sbatch`, job 32986701)

Stages `tools fmd runs` of `rebuild_hprc2.sh`: download `human472.fmr.gz` from Zenodo 14854401, check its
md5, convert it with the pinned ropebwt3 (`build -i -do`), and run fmd2rlbwt.

All three equal the frozen files (`cheap-result.txt`): `ropebwt3/human472/index.fmd` (sha256 10e89eb1...),
`precursors/human472/bwt.heads` (7814fdcd...) and `precursors/human472/bwt.len` (aa5d481e...). The `.fmr.gz`
passed the md5 Li's record lists.

| Stage | Wall | Peak RSS |
|---|---|---|
| tools (clone and build fmd2rlbwt and TeraLCP) | 1 min 13 s | 0.3 GB |
| fmd (8.0 GiB download, gunzip, `ropebwt3 build -t 8 -i -do`) | 27 min 40 s | 40 GB |
| runs (fmd2rlbwt) | 2 min 12 s | 0.08 GB |

## Not checked

The HPRC2 `min` and `phi` (stage `lcp`, about 13 h and 171 GiB) and the s8 and s6 SSAs (stage `ssa`, about
10 h each and 413 GiB of disk) were not rebuilt. The medium tier shows that the same TeraLCP build
reproduces the smaller deposited precursors.
