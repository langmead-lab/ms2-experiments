# Fig rep's counting runs outside the freeze

Fig rep plots 16 batches. Twelve come from the freeze's A1 cells (`fig_rep-a2` and `fig_rep-b2` in
`freeze-2885ee8-00e665de`): mtb152 at 0 to 0.5% errors and the chr22 indexes of 5 to 100
haplotypes. The other four were added after the freeze, on Rockfish, and their records are here:

| Run | Host | Started (UTC) | Reads |
|---|---|---|---|
| `mtb152-A2-e0.010` | devlangmead1 | 2026-11-01 14:24 | 20 Mbases from the 8 held-out mtb152 genomes, 1% errors (`inputs/mtb152-A2-e0.010/manifest.json`) |
| `mtb152-A2-e0.020` | devlangmead1 | 2026-11-01 14:24 | as above, 2% errors |
| `chr22-B-h200` | sr06 (Slurm) | 2026-11-01 16:29 | the freeze's simulated chr22 PacBio HiFi reads (sha256 `8e5897988e47...`, the same file as the freeze's chr22-B cells) |
| `chr22-B-h464` | sr06 (Slurm) | 2026-11-01 17:09 | as above |

Each `runs/<run>/run.json` records the command lines, binary, index paths, reads and their sha256.
`summary.json` is the per-batch summary that `fig_rep_tikz.py` reads, in the same format
(`fig-rep-summary-v1`) as the freeze cells' `rep/summary.json`.

## Build

The binary is `ms2-replog`, the instrumented counting build of ms2 at 2885ee86 (`paper-v1`), built on
Rockfish (sha256 `d70dc5f37b75...`; the freeze's AWS build of the same commit is `2b123fdd3470...`).
The mtb152 indexes were rebuilt on Rockfish from the same precursors as the freeze's; they have the
freeze indexes' sizes but not their sha256. The chr22 indexes of 200 and 464 haplotypes are new.

## Cross-check

`runs/mtb152-A2-e0.0050.check` reran the freeze's 0.5% batch with this binary and the Rockfish
mtb152 indexes, on the freeze's reads (sha256 `4351e2c79930...`). Its summary equals the freeze
cell's in every field but the timings, so the counts do not depend on which build or index copy ran.
The figures driver leaves this run out of Fig rep.

## Not here

The `.replog` files that the violins are drawn from stay on Rockfish, under
`/scratch4/blangme2/ms-bakery/figrep-axes/runs/`. The box statistics and means that Fig rep's
ticks, boxes and caption use are all in `summary.json`.

## Code

`code/` holds the scripts as they ran: `fig_rep.py` (runs the three modes and writes `run.json`
and `summary.json`), `prepare_more.py` and `plan_more.py` (make the 1% and 2% read sets and
plan the new chr22 indexes; `plan-more.json` is the plan), `build_cell.sbatch` and `run_b.sbatch`
(the Slurm jobs for the chr22 builds and runs), and `fig_rep_violin.py` and `rep_zero.py` (the
violin densities and the zero-operation fractions).
