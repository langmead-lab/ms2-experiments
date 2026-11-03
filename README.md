This repo contains the code and records behind the experiments and measurements for the paper
titled "Fast Pangenome Matching Statistics in Compressed Space with Move Structures" by Nathaniel K.
Brown, Ahsan Sanaullah, Pramesh Shakya, Shaojie Zhang, and Ben Langmead. The `ms2` software itself
is at [langmead-lab/ms2](https://github.com/langmead-lab/ms2). The paper used the version at commit
[2885ee86](https://github.com/langmead-lab/ms2/commit/2885ee86270b4187fc45dcc69efd7de48a4aac03),
which is tagged `paper-v1`.

| Path | Content |
|---|---|
| `ms2exp/` | Benchmark workflow Python package; uses standard library only, except figure rendering |
| `launchers/rockfish/`, `launchers/aws/` | Two "launchers" we used: one using Slurm on the JHU Rockfish cluster, and another for AWS spot instances (CDK app and worker; copy `launchers/aws/config.example.json` to `config.json` first) |
| `results/` | Small recorded measurements from the two freeze runs; `results/README.md` maps them to paper elements |
| `figures/` | Output from `python -m ms2exp.figures`: summaries, drafts, TikZ data and the paper's figure PDFs |
| `tools.lock.json`, `TOOLS.md` | The tools pinned to specific versions, and how each was built |
| `patches/` | The fixes to PHONI and Augmented PHONI that make their matching statistics correct (`TOOLS.md`) |
| `toy/` | A small "toy" end-to-end run of the whole workflow, run directly on one machine (no Slurm) |
| `data/` | Input download and preparation scripts; `data/locations/` maps logical input names to Rockfish paths (the Langmead lab's own copies) |
| `docs/` | Design records of the freeze and thread-scaling runs. Code docstrings and the committed result records cite them by bare file name: `BENCH_FREEZE_DESIGN.md` (decisions D1 to D13) and `THREAD_SCALING_PLAN.md` |
| `tests/` | Unit tests (`python -m pytest tests`) |

## Workflow

The workflow has three parts. The **manifest** (`ms2exp.manifest`) lists jobs and cells and
"pinned" inputs, tool binaries included, by sha256, giving each a logical name. A **cell** is one tool
configuration run on one read set.  All measurements run the combination three times with a "cold"
cache, with run 1's output checked against true statistics. The **runner** (`ms2exp.run.runner`) runs
one job on either platform: it stages and verifies inputs, runs cells in a random order, checks them,
and writes per-cell records once it finishes.  An interrupted job can resume where it stopped. One
**launcher** per platform starts runners: `launchers/rockfish` for staging, dry runs and counting
runs (i.e. runs that produce needed outputs but aren't important to time/measure), and `launchers/aws`
for the timed runs. `ms2exp.collect` then gathers the cells into one summary per cell group and checks
them.

| Module | Role |
|---|---|
| `ms2exp.lib` | Manifests, hashing, local and S3 stores |
| `ms2exp.catalog` | Inputs catalog: `hash`, `link`, `publish` |
| `ms2exp.manifest` | `--matrix freeze` or `--matrix toy` |
| `ms2exp.seed_rerun` | A rerun manifest under a changed protocol (the THP rerun) |
| `ms2exp.run.runner`, `.timing`, `.common` | The runner, phase timers and counting-build stats, command execution |
| `ms2exp.validate.ms`, `.offsets`, `.rope`, `.legacy` | Truth checks for ms2 and ropebwt3 MS, locate output, ropebwt3 SMEMs, MONI, PHONI and Augmented PHONI |
| `ms2exp.build.tools`, `.job`, `.d13_prep` | Tool builds, AWS build jobs, D13 read-set preparation |
| `ms2exp.sim.chr22_sim` | Simulated read sets and their truth |
| `ms2exp.collect` | Group summaries and cross-tool checks |
| `ms2exp.figures` | Figure and table generators, and the driver below |

## Figures

```bash
python -m venv .venv && .venv/bin/pip install -e '.[figures,test]'
.venv/bin/python -m ms2exp.figures --results results --out figures
```

This rebuilds `figures/` from `results/freeze-2885ee8-00e665de`. It needs latexmk for the PDFs
(`--no-pdf` skips them). The PDFs in `figures/pdf/` carry the names the paper includes:

| Paper figure | File | TikZ source |
|---|---|---|
| Fig ops | `ops.pdf` | `fig-ops-row.tex` |
| Fig ms22 | `ms22_100.pdf` | `fig-ms22-100.tex` |
| Fig ms-oracle | `ms_oracle.pdf` | `fig-ms-oracle.tex` |
| Fig locate100 | `locate100.pdf` | `fig-l100.tex` |
| Fig ms22-10 (appendix) | `ms22_10.pdf` | `fig-ms22-10.tex` |
| Fig locate1, locate10 (appendix) | `locate1.pdf`, `locate10.pdf` | `fig-l1.tex`, `fig-l10.tex` |
| Fig ms-ftab (appendix) | `ms_ftab.pdf` | `fig-ms-ftab.tex` |
| Fig locate-ftab (appendix) | `locate_ftab.pdf` | `fig-lftab.tex` |
| Fig threads (appendix) | `threads.pdf` | `fig-threads.tex`, from `results/threads-2885ee8-da1041e8` |

Every number plotted is in `figures/draft/plotted-points.csv`, with its median query and load
times. The load-time table is in `figures/draft/CAPTIONS.md`. Fig rep is not regenerated from these
results (`results/README.md` says why).

## Small "toy" run

A good test to see if this repo is basically working in your environment is to run:
`toy/toy_prep.sh TOY_DIR MS2_REPO MS2_COMMIT [--submodule-repo NAME URL ...]`.  This builds the tools
and small mtb152 indexes at the given commit of the `ms2` repo, creates small read sets, and writes a
small manifest whose cells cover every path of the runner. `toy/toy_run.sh TOY_DIR` runs it and
collects.  Both currently assume you are running on the JHU Rockfish cluster, so some updates will
certainly be needed to get it working on yours. While ms2 is private, pass git bundles of the repository
and both submodules, made by `toy/make_source_bundles.sh`.

## Tools

`tools.lock.json` pins ms2 and its submodules, ropebwt3, MONI, PHONI and Augmented PHONI, and the
runtime libraries shipped with them; `TOOLS.md` says how each binary was built and checked.

## Acknowledgements

`ms2` builds on Orbit, the move-structure library by Nathaniel Brown, and TeraTools, the move-structure building library and its TeraMS component for computing matching statistics. Its psi and phi query engine follows the design and code of TeraTools (TeraMS and TeraLCP) by Ahsan Sanaullah, Shaojie Zhang and the S. Zhang Lab at UCF, and it builds from TeraLCP's output. fmd2rlbwt uses ropebwt3's FMD reader by Heng Li.

## License

Software in this repo is under the MIT license (`LICENSE`). The patches in `patches/` are under their upstream
projects' licenses.
