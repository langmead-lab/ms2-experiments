# Inputs

This is the entry point for every input of the freeze manifests (`results/freeze-2885ee8-00e665de/manifest.json`, key `inputs`; the sha256 of each file is there). Each group below gives its public origin, the steps that derived it, and where to get it. "Deposit 1, 2, 3" are the three Zenodo records: [10.5281/zenodo.23964218](https://doi.org/10.5281/zenodo.23964218), the read sets with truth; [10.5281/zenodo.23964220](https://doi.org/10.5281/zenodo.23964220), the chr22 and mtb152 precursors with the chr22 and mtb152 ropebwt3 FMDs; [10.5281/zenodo.23964222](https://doi.org/10.5281/zenodo.23964222), the MONI, PHONI and AugPHONI chr22 indexes. Each record's file names are the manifest paths with `/` replaced by `__`, and its `SHA256SUMS` matches `inputs.tsv`; `python3 -m ms2exp.fetch` downloads and checks them. "Rebuild" means the file is not deposited and has to be made again.

Terms: an **FMD** is ropebwt3's FM-index of both strands; the **RLBWT precursors** are the run heads and lengths of the Burrows-Wheeler transform (`bwt.heads`, `bwt.len`), the per-run LCP minima (`min`) and the phi samples (`phi`) that ms2 builds its indexes from; **truth** (`truth.u32`) is the matching statistic at every query position; an **SSA** is ropebwt3's sampled suffix array.

| Directory | Contents |
|---|---|
| `download/` | Source downloads: `download.py` (ENA runs in `manifests/ena-files.json`), `prepare_hg002.py` (the Nanopore sample), `download_skin.py` (the skin run, or all skin runs of the study) |
| `prep/` | Drivers for the preparation steps, and the reference and index builds |
| `lists/` | The 472 HPRC2 haplotypes, the HPRC2 assembly order and the chr22 panel |
| `evidence/` | Records of how each frozen file was made on Rockfish (`evidence/README.md`) |
| `locations/` | Rockfish paths of the frozen inputs (for the Langmead lab) |

The preparation modules are in `ms2exp/prep/`, the simulators in `ms2exp/sim/chr22_sim.py`, and the read-set steps of the freeze in `ms2exp/build/d13_prep.py`. Some steps still name Rockfish paths as their defaults; each script's usage line says how to point it elsewhere. Tool versions and sha256s are in `tools.lock.json` (`prep_tools` for the preparation tools).

## References and indexes

| Group | Manifest paths | Public origin | Derivation | Obtain |
|---|---|---|---|---|
| HPRC2 haplotypes | (content of every HPRC2 index) | `human472.agc` in Zenodo 14854401 (Li): the February 2025 HPRC release 2 data (r2 v0.6, 464 haplotypes) plus CHM13, GRCh38, CN1, YAO and KSA001 | `lists/human472-haplotypes.txt` (472 AGC sets with PanSN prefixes); `python -m ms2exp.prep.lists haplotypes --len human472.fmd.len.gz` checks it against Li's index | This repository |
| ropebwt3 HPRC2 FMD, lengths, s11 SSA | `ropebwt3/human472/index.fmd`, `.len.gz`, `.s11.ssa` | Zenodo 14854401 (`human472.fmr.gz`, `human472.fmd.len.gz`, `human472.fmd.ssa.gz`) | `ropebwt3 build -i human472.fmr -do human472.fmd`; the length file unchanged; the SSA decompressed | Li's record |
| ropebwt3 HPRC2 s8, s6 SSAs | `ropebwt3/human472/index.fmd.s8.ssa`, `.s6.ssa` | derived | `ropebwt3 ssa -s8`, `-s6` on the FMD | Rebuild: `prep/rebuild_hprc2.sh WORK fmd ssa` (about 10 h each) |
| HPRC2 precursors | `precursors/human472/*` | derived from the FMD | `fmd2rlbwt` (heads, lengths), `TeraLCP -f fmd -ominima -ophi` (TeraTools 414f809, tag `ms2-precursors`) | Rebuild: `prep/rebuild_hprc2.sh WORK tools fmd runs lcp` (about 13 h and 171 GiB for TeraLCP) |
| chr22 references | (content of every chr22 index) | `human472.agc`; chr22 assignments by RagTag 2.1.0 against CHM13 (`evidence/ragtag/`) | `prep/chr22_panel.py extract N` from `lists/chr22-panel.tsv`: the first N haplotypes of `lists/hprc2-donor-order.txt` without CHM13 and GRCh38 | Rebuild (minutes) |
| chr22 ropebwt3 FMDs | `ropebwt3/chr22-h10/index.fmd`, `chr22-h100/...` | derived | `prep/legacy_indexes.py --tools ropebwt3` (`ropebwt3 build -t2 -d`) | Deposit 2 |
| chr22-h10, h100 precursors | `precursors/chr22-h{10,100}/*` | derived from the chr22 FMDs | `prep/chr22_precursors.sh` (fmd2rlbwt, TeraLCP bce9ea7 minima renumbered by fmd2rlbwt, TeraLCP 414f809 phi) | Deposit 2 |
| A1 ladder precursors | `precursors/chr22-h{1,2,5,20,50}/*` | derived from the chr22 references | `prep/ladder_precursors.sh` (`build_rlbwt.cpp` suffix sort of the %-joined text, TeraLCP af42783 `-f rlbwt` minima) | Deposit 2 |
| mtb152 split, FMD and precursors | `ropebwt3/mtb152/index.fmd`, `precursors/mtb152/*` | Li's ropebwt3 index of 152 M. tuberculosis genomes (`mtb152.tar.gz`) | `prep/mtb152_split.sh`, reconstructed from Nate Brown's log: the first 8 genomes held out, the other 144 rebuilt with `ropebwt3 build -do`, then fmd2rlbwt and TeraLCP 414f809 | Deposit 2 |
| MONI, PHONI, AugPHONI indexes | `moni/*`, `phoni/*`, `augphoni/*` | derived from the chr22 references | `prep/legacy_indexes.py --tools moni phoni aug_phoni` (tool builds in `TOOLS.md`) | Deposit 3 |
| ms2 indexes | `ms2/*` | derived from the precursors | built by the freeze manifest's build jobs | Rebuild by the freeze |

## Read sets

Every 60 Mbase set holds whole reads made only of A, C, G and T, taken in sampled order up to at most 60,000,000 bases. The locate sets are whole-read prefixes of them. `prep/run_d13.sh` repeats the last steps (caps, samples, prefixes, the chr22 simulations and the skin set) in the order they ran.

| Group | Manifest paths | Public origin | Derivation | Obtain |
|---|---|---|---|---|
| CHM1 short reads | `reads/human472/human_short-*` | ENA SRR3099549 | `download/download.py`; `ms2exp.prep.pool human_short` (SeqKit 2.10.1 `sample -s 20261004 -p 0.01737` on both mate files, then `-s 20261005 -p 0.2` for the 200 Mbase core) and `ms2exp.prep.truth_pool`; `d13_prep.py short-60mb`; `cap --limit 6000000` | Deposit 1 |
| PacBio HiFi HG00099 | `reads/human472/human_pacbio-*` | ENA SRR26545347 (HiFi, DeepConsensus v1.2) | `download.py`; `ms2exp.prep.pool human_pacbio` (`-s 20261004 -p 0.006004`, then the core) and `truth_pool`; `ms2exp.prep.acgt_workload core`; `d13_prep.py cap`; `cap --limit 20000000` | Deposit 1 |
| Nanopore HG002 | `reads/human472/human_ont_modern-*` | ONT open data gm24385_2023.12 (`all_pass.vhg002v1.bam`) | `download/prepare_hg002.py` (samtools `--subsample 0.008 --subsample-seed 20261003`); `ms2exp.prep.pool human_ont_modern` (the whole sample, then the core) and `truth_pool`; `acgt_workload core`; `cap`; `cap --limit 10000000` | Deposit 1 |
| Skin metagenome | `reads/human472/guccione_skin-*` | ENA ERR14066769 (PRJEB83637, Guccione et al. 2025) | `download/download_skin.py`; chosen among the study's skin runs by `ms2exp.prep.skin_screen` (`evidence/skin/`); `d13_prep.py skin-60mb`; `cap --limit 500000` | Deposit 1 |
| APR008 chr22 long arm | (source of the chr22 simulated reads) | GenBank GCA_050494745.1, scaffold JBHMUJ010000002.1 (34,044,963 bases, no N) | `prep/apr008_chr22_map.sh` (minimap2 2.28 against CHM13 v2.0 chr22; `evidence/rockfish/results/freeze-d13/apr008-chr22/`); `d13_prep.py apr008-chr22` | NCBI |
| chr22 simulated reads | `reads/chr22-h{10,100}/*`, `reads/chr22-sim/*` | derived | `d13_prep.py simulate --set short|pacbio|ont` (the `ms2exp/sim/chr22_sim.py` recipe and seeds, truth on each chr22 FMD) | Deposit 1 |
| Genome pieces | `reads/human472/human472-apr008-hap1-p100k-s600-v1`, `reads/mtb152/mtb152-heldout16-p100k-s600-v1` | APR008 hap1 (all 79 scaffolds); the 8 held-out mtb152 genomes with their reverse complements | `prep/fig_ops_sets.sh` (`ms2exp.prep.acgt_split`, `ms2exp.prep.pieces workload`, `truth`, `pieces --length 100000`); `d13_prep.py sample --seed 2026103001` (HPRC2, 600 of 30,424 pieces) and `2026103002` (mtb152, 600 of 718) | Deposit 1 |
| Fig ops PacBio reads | `reads/human472/human472-apr008-pacbio-60mb-v2`, `reads/mtb152/mtb152-pacbio-60mb-v2` | the same held-out genomes | `prep/fig_ops_sets.sh` (`chr22_sim simulate --set pacbio --source ... --seed 2026102612` for HPRC2, `2026102611` for mtb152, then `pieces workload` and `truth`); `d13_prep.py cap` | Deposit 1 |
| A1 error series | `reads/mtb152/mtb152-A2-*-v1` | the 8 held-out mtb152 genomes | `prep/series_a2.sh` (`ms2exp.figures.fig_rep simulate`, seeds 2026102801 to 2026102810) | Deposit 1 |

The pool steps (`pool`, `truth_pool`, `acgt_workload`) ran before the freeze and are kept as provenance of the frozen pools; they are not part of `run_d13.sh`, which starts from those pools. In the CHM1 set, both mates of a pair were sampled together by the first SeqKit step, but the nested 20% core sample drew reads independently, so the set holds some single mates. Each read was a separate query.

## Truth

Each `truth.u32` holds the matching statistic at every query position. It came from ropebwt3 (binary sha256 cfad8475, built from commit b408339) run as `ropebwt3 mem -l1` on the FMD of the index being queried: the human472 FMD for the real reads (on the 200 Mbase cores, then sliced to each set), the skin set and the HPRC2 genome pieces; the chr22-h10 and chr22-h100 FMDs for the simulated chr22 reads; and the mtb152 95% FMD for the mtb152 sets. `ms2exp.validate.ms.check` converted the super-maximal exact matches to matching statistics, and every value was checked to lie between 1 and the distance to the read end. For the genome pieces, the truth of the whole held-out genome was clamped at each piece end, min(MS[i], end - i). The A1 replay runs (the error series and the chr22 ladder) have no truth; they require psi, phi and dual to produce identical output.

ropebwt3 is also one of the benchmarked tools. Agreement of the ms2 outputs, and of the corrected PHONI and AugPHONI outputs, with this truth at every position is the cross-check; for the real-read cores, the Orbit and TeraMS runs of the real-hprc2 campaign also matched it at every position (`evidence/rockfish/results/campaign-20261012/inputs/*/explore-acgt-60mb-v1/manifest.json`).
