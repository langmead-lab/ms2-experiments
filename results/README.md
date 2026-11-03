# Results of the freeze benchmarks

Two manifests produced every timed number in the paper apart from the thread-scaling figure, whose
manifest is described in "Thread scaling" below. Both ran ms2 at 2885ee86 (`paper-v1`) on
AWS spot instances (i4i.metal), one thread per command, under `numactl` on socket 0, three cold
runs per cell with the median reported, and run 1 of every cell checked against truth.

| Manifest | Run | Role |
|---|---|---|
| `freeze-2885ee8-3c8c8967` | Oct 30 to 31, 2026 | Transparent huge pages (THP) off: the protocol did not set THP and the host records do not include the setting, which the runner records from the rerun on. The THP attribution compares the 230 cells timed in both runs. Of the 279 cells this run recorded, the other 49 (count and Fig rep cells) were copied into the rerun rather than rerun. The 15 ms2-psi locate cells did not finish in this run |
| `freeze-2885ee8-00e665de` | Oct 31 to Nov 1, 2026 | The same jobs, cells and inputs with THP set to `always` (`seed_rerun`); every figure and table of the paper |

Each has 83 jobs and 294 cells. The rerun copied the first run's build records, oracle traces and
untimed jobs (`seeded-from.json`), so only the timed cells ran again.

## What is here

Per manifest directory: `manifest.json` (the manifest, with its inputs catalog pinning every input
by sha256), and per job `JOB/cells/CELL/` with `summary.json`, the run records (`run1.json` ...),
each command's `*.command.json` and stderr, and the `CELL.done` marker holding the summary's
sha256; `JOB/attempts/` and `JOB/job.done`; build jobs' `outputs.json` (the sha256 and size of
every built index and binary); `boxes/` (the worker log of every AWS box); `seeded-from.json` in the
rerun. About 15 MB and 3,300 files per manifest.

Not here, but recorded by key, size and sha256 in each directory's `s3-excluded.tsv`: the locate
outputs in canonical form (`canonical.gz`), the oracle traces (`.oracle`) and the repositioning logs
of Fig rep (`.replog`), about 24 GB in all. They stay in the lab's freeze bucket, `results/` of the
S3 root that the launcher config names (`bench_s3_root`; the as-run configs are in `results/code/`). The sha256 of each canonical
file's uncompressed content is also in its cell's `run1.json` (`check.canonical.sha256`).

## From results to figures

`python -m ms2exp.figures` (see the top-level README) runs `ms2exp.collect --platform aws
--no-compare` on `freeze-2885ee8-00e665de` and writes one summary per cell group to
`figures/summaries/`. Groups and figures:

| Cell groups | Paper figure (`figs/` name) |
|---|---|
| `fig_ops-*` and their `-counts` groups | Fig ops (`ops`) |
| `chr22-h100-ms-sim*` | Fig ms22 (`ms22_100`) |
| `chr22-h10-ms-sim*` | Fig ms22-10, appendix (`ms22_10`) |
| `*-ms-acgt60mb` | Fig ms-oracle (`ms_oracle`); the short-read group also Fig ms-ftab (`ms_ftab`) and Table load |
| `*-locate-*` | Fig locate100, locate1, locate10 (`locate100`, `locate1`, `locate10`); the short-read group also Fig locate-ftab (`locate_ftab`) |
| `fig_rep-a2`, `fig_rep-b2` | Fig rep's A1 runs (see below) |

`--no-compare` skips the locate cross-tool comparison, which needs the canonical files. The
comparison ran when the rerun was collected with the canonical files present: all 135 pairs of an
ms2 locate cell and a ropebwt3 cell of the same read set and cap agreed. The D9 check (dual on the full index with sequence names against
the plain full index) passes in every group: query times within 0.6%.

Fig rep is the exception. The paper's current Fig rep data (`fig_rep_tikz.py --design 5`) include
rungs the freeze never ran (1% and 2% error, 200 and 464 haplotypes), so they come from runs
outside these manifests. `ms2exp.figures` writes the freeze's own A1 series (`--design 4`, box
statistics only, since the violins need the `.replog` files).

## The code that ran

The AWS boxes of the two freeze manifests ran the freeze code from a tarball, `code/freeze-code.tar.gz`
under the S3 root, in the ms-bakery layout from which `ms2exp` was later extracted (commit 0ff1ebb).
Each box downloaded the tarball when it booted. Every version is in `results/code/`, with its sha256
in `results/code/SHA256SUMS`:

| File | Ran | In S3 (UTC) | sha256 |
|---|---|---|---|
| `freeze-code.before-fig3count-fix.tar.gz` | `freeze-2885ee8-3c8c8967` until the fig3count fix | until 2026-10-31 12:03, then kept under this name | 7de51427d97f3a701a4d914ba88de3455cf8d18c145f6f0bf7154a8cae31583f |
| `freeze-code.before-thp-rerun.tar.gz` | `freeze-2885ee8-3c8c8967` after the fix | 2026-10-31 12:03 to 21:49, then kept under this name | 18f16c6c15652a118a6f1c1f395ad08ef763844cf6e6ec43154fd05d6880390f |
| `freeze-code.tar.gz` | `freeze-2885ee8-00e665de` | from 2026-10-31 21:49 | 300e2338a0a63d689359e4954478225f327e01df276dead353f8365b893132d6 |
| `threads-code.tar.gz` | `threads-2885ee8-da1041e8` | 2026-11-02 22:32, as `threads/code/freeze-code.tar.gz` | 1ddfd4ee5806e3277bdb0a798472652b0517124df1ff8c290c39e79bb04c873e |

### The fig3count fix

The operation counts of Fig ops come from the counting build `ms2-stats`, in the four `fig3count-*`
jobs. In the first version of the freeze code, two bugs stopped some of their cells from finishing.
Neither could produce a wrong number.

- The untimed preparation of a cell recorded an oracle trace only when the cell's tool was `ms2`, so
  `ms2-stats` cells in oracle mode started without a trace and exited with an error. The box logs
  record 168 such failures between 02:28 and 11:58 UTC on October 31.
- The parser of the counting build's `stats:` lines took a line named for a mode (`stats: dual ...`),
  which repeats `psi_steps/base`, as a second operations line, and refused a dual-mode cell with
  "expected one stats line".

The fix records the trace for every ms2 tool, or fetches the copy that another job already stored,
and files the mode-named lines under the record's `summary_lines`. Two unrelated changes went with
it: the worker's check for stale job claims, and the `libstdc++-static` package that the boxes
install for the static ms2 build. The first fixed oracle cell ran at 12:08 UTC on a box that had
booted at 00:33, so the fixed code reached boxes that were already running.

Of the 32 count cells, 29 finished after the fix, between 12:08 and 13:01 UTC on October 31. The
other three (`mtb152-heldout16` psi and phi ordinary, `mtb152-pacbio` psi ordinary) finished between
02:29 and 02:32 UTC with the first version. Parsing every committed `run1.stderr` of the count cells
again with both versions of the parser gives the committed operation counts (`ops_per_base`) for all
64 count records of the two manifests. On the early phi record the fixed parser adds only a
`summary_lines` entry, which no figure reads. The rerun copied the four count jobs from the first run
unchanged (`seeded-from.json`).

### The THP rerun

`freeze-code.tar.gz` differs from its predecessor only in transparent huge pages and in the new
`seed_rerun.py`: the protocol sets THP to `always`, the box sets it at boot, every timed job checks it
before it starts, and the host record includes it.

### Equivalence with `ms2exp`

The toy run (`toy/`) checks that `ms2exp` behaves as the snapshot did. On 2026-11-01 both ran the
toy on devlangmead1 (the lab's development server on Rockfish) at ms2 2885ee86 from the same source bundles, and `toy/compare.py` found the
same group summaries, collect report, cell records and done markers, apart from timings, memory,
hosts and paths. The only other differences were the sha256s of the separately built toy indexes,
whose headers record their build path and time. The `ms2exp` runner and collect, run on the
snapshot side's own manifest and inputs, gave identical records. `toy/expected-v1-snapshot/` holds the
snapshot side's summaries; `toy/expected/` holds those of the current toy, on the v2 read sets. The toy covers the mtb152 and chr22-h10 cells; the HPRC2 cells run the
same code at larger sizes.

After that check, commit 49d08f0 added thread scaling to `ms2exp/run/runner.py`,
`ms2exp/run/common.py`, `manifest.py`, `collect.py`, `catalog.py` and `d13_prep.py`. A cell without
a `threads` field, which includes every freeze cell, keeps the freeze command line, NUMA binding and
environment (`tests/test_threads.py`). `collect.py` adds a `cross_thread` key to `report.json`, which
is empty for the freeze manifests. `paper-v1` includes these changes.

## Thread scaling

One more manifest produced the thread-scaling figure (Fig threads, `fig:threads` in the appendix).
It ran the freeze binaries and indexes (ms2 2885ee86, ropebwt3 b408339) on AWS i4i.metal spot
instances, one job per configuration and read set, each on a box of its own. Its protocol is the
freeze protocol with THP `always`, plus a thread count per cell: T = 1 to 32 bound to the first T
physical cores of socket 0, and T = 64 on both hardware threads of those 32 cores, with memory on
node 0. Every run records the sha256 of its output, and `collect.py` checks that each configuration
gives one output at every thread count.

| Manifest | Run | Role |
|---|---|---|
| `threads-2885ee8-da1041e8` | 2026-11-02 22:48 to 2026-11-03 02:51 UTC; 11 jobs, 71 cells | Fig threads: ms2 min, psi, phi and dual (2 MiB blocks) and ropebwt3 `-l1` on the HPRC2 short-read and ONT 1 Gbase sets, at 1, 2, 4, 8, 16, 32 and 64 threads, plus one 64-thread check cell on a 60 Mb set |
| `threads-pilot-2885ee8-5a86a009` | 2026-11-02; 6 jobs, 21 cells | The pilot, at 1 and 32 threads on the Fig 5 60 Mb sets (ms2 on one AWS box, ropebwt3 on Rockfish); in no figure |

All 71 cells of `threads-2885ee8-da1041e8` passed their truth checks, and each of its 11 groups of a
configuration and read set gave one output sha256 at every thread count (`report.json` there).
`figures/tikz/fig_threads_data.py` reads the cell summaries of the 1 Gbase cells and writes
`figures/tikz/data/fig-threads/*.dat`, which `figures/tikz/fig-threads.tex` plots; the caption is
`figures/tex/fig-threads-caption.tex`. The boxes ran `results/code/threads-code.tar.gz`, which holds
the same `ms2exp` and launcher files as `paper-v1`.
