# Preparation evidence

Small records copied from the lab's Rockfish cluster so that each input's derivation can be read without lab access. Nothing here is needed to run the benchmark. `data/README.md` says which record backs which input.

## `rockfish/`

Files under `rockfish/` keep their path relative to `/scratch16/blangme2/ms-bakery` on Rockfish.

| Path | What it records |
|---|---|
| `results/real-hprc2-1d580b6/inputs/COHORT/complete.json` | The SeqKit sampling of each real read pool (commands, seeds, fractions, counts and sha256 of every set), written by `ms2exp/prep/pool.py`'s original. |
| `results/real-hprc2-1d580b6/truth/COHORT/core.json` | The ropebwt3 `mem -l1` truth run on each 200 Mbase core pool, with its whole-output check. |
| `results/campaign-20261012/inputs/COHORT/explore-acgt-60mb-v1/manifest.json` | The ACGT-only whole-read prefix of the PacBio and Nanopore cores (`ms2exp/prep/acgt_workload.py`), including the other tools that matched the core truth. |
| `results/campaign-20261012/inputs/heldout/apr-APR008-hap1/` | The APR008 haplotype 1 download (URLs, publisher checksums), the evidence that APR008 is not in HPRC2, and the original `acgt_split.py`. |
| `results/campaign-20261012/inputs/fig_ops/`, `fig_ops2/` | Genome query workloads, their 100 kbase piece sets and the simulated PacBio sets from the held-out genomes (`ms2exp/prep/pieces.py`, `ms2exp/sim/chr22_sim.py`), with truth records. |
| `results/campaign-20261012/fig_rep/inputs/mtb152-A2-*/manifest.json` | The A1 mtb152 error series (`data/prep/series_a2.sh`). |
| `results/freeze-d13/` | The D13 read sets as frozen: the original Slurm drivers (`prep-*.sbatch`, reproduced by `data/prep/run_d13.sh`), the per-set manifests and truth records, and `apr008-chr22/`. |
| `results/freeze-d13/apr008-chr22/` | The APR008 chr22 identification: `map.sbatch` (reproduced by `data/prep/apr008_chr22_map.sh`), `minimap2.stderr`, `per-scaffold.tsv` and `long-arm.tsv`, and the sha256 of the 96 MB PAF. |
| `data/ladder/chr22-hN/manifest.json` | Each chr22 reference: its AGC components, extracted-FASTA sha256s and ACGT segment coordinates (`data/prep/chr22_panel.py`). |
| `freeze-catalog/chr22_phi.sbatch` | The phi samples of chr22-h10 and h100 (reproduced in `data/prep/chr22_precursors.sh`). |

`per-scaffold.tsv` has no header. Its columns are scaffold, scaffold length, query bases in alignments of mapping quality 20 or more, and the first and last CHM13 chr22 coordinate those alignments cover. It sums every alignment, so scaffolds of other chromosomes also show tens of megabases. `long-arm.tsv` keeps only alignment blocks of at least 500 kbases (PAF column 11) at mapping quality 20 or more. It shows that JBHMUJ010000002.1 aligned 33,878,807 bases in two blocks covering CHM13 chr22 17,473,817 to 51,324,621, and that every other scaffold's blocks of that size lay below 17.4 Mbases.

## Other records

| Path | What it records |
|---|---|
| `fmdin/stageA.sh` | The run dump and LCP minima of chr22-h10 and h100 (`/vast/blangme2/hprc2-ms/dev/fmdin/stageA.sh` on Rockfish), reproduced in `data/prep/chr22_precursors.sh`. |
| `ragtag/` | How the chromosome assignments behind `data/lists/chr22-panel.tsv` were made: the lab's RagTag job template and driver (`/data/blangme2/fasta/hprc-v2/ragtag/`), and the first lines of one job's log, which give RagTag v2.1.0 and its command. |
| `mtb152/setup.20261001_132202.log` | Nate Brown's log of the mtb152 95/5 split (ropebwt3 version lines removed), the only record of the commands that `data/prep/mtb152_split.sh` reconstructs. |
| `skin/guccione-screen-v1.json` | The screen that chose skin run ERR14066769 (`ms2exp/prep/skin_screen.py`). |
