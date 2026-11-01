# Results of the freeze benchmarks

Two manifests produced every timed number in the paper. Both ran ms2 at 2885ee86 (`paper-v1`) on
AWS spot instances (i4i.metal), one thread per command, under `numactl` on socket 0, three cold
runs per cell with the median reported, and run 1 of every cell checked against truth.

| Manifest | Run | Role |
|---|---|---|
| `freeze-2885ee8-3c8c8967` | Oct 30 to 31, 2026 | Transparent huge pages (THP) off: the protocol did not set THP and the host records do not include the setting, which the runner records from the rerun on. The THP attribution compares its 279 comparable cells with the rerun |
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
S3 root that `launchers/aws/config.json` names (`bench_s3_root`). The sha256 of each canonical
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

The AWS boxes ran the freeze code from a tarball, `code/freeze-code.tar.gz` under the S3 root, in
the ms-bakery layout that `ms2exp` was extracted from. Its versions:

| File | Uploaded | sha256 |
|---|---|---|
| `freeze-code.before-fig3count-fix.tar.gz` | 2026-10-31 12:03 UTC | 7de51427d97f3a701a4d914ba88de3455cf8d18c145f6f0bf7154a8cae31583f |
| `freeze-code.before-thp-rerun.tar.gz` | 2026-10-31 21:49 UTC | 18f16c6c15652a118a6f1c1f395ad08ef763844cf6e6ec43154fd05d6880390f |
| `freeze-code.tar.gz` (ran `freeze-2885ee8-00e665de`) | 2026-10-31 21:49 UTC | 300e2338a0a63d689359e4954478225f327e01df276dead353f8365b893132d6 |

The toy run (`toy/`) checks that `ms2exp` behaves as the snapshot did. On 2026-11-01 both ran the
toy on devlangmead1 at ms2 2885ee86 from the same source bundles, and `toy/compare.py` found the
same group summaries, collect report, cell records and done markers, apart from timings, memory,
hosts and paths. The only other differences were the sha256s of the separately built toy indexes,
whose headers record their build path and time. The `ms2exp` runner and collect, run on the
snapshot side's own manifest and inputs, gave identical records. `toy/expected/` holds the
snapshot side's summaries.
