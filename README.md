# ms2-experiments

The code and records behind every timed number, figure and table of the ms2 matching-statistics
paper. ms2 itself is [langmead-lab/ms2](https://github.com/langmead-lab/ms2); the paper used commit
2885ee86, tagged `paper-v1`.

| Path | What it holds |
|---|---|
| `ms2exp/` | The benchmark workflow as one Python package (standard library only, except the figures) |
| `launchers/rockfish/`, `launchers/aws/` | The two platforms' launchers: Slurm on Rockfish, and an AWS spot fleet (CDK app and worker) |
| `results/` | The small records of the two freeze runs; `results/README.md` maps them to the paper |
| `figures/` | Everything `python -m ms2exp.figures` writes: summaries, drafts, TikZ data and the paper's figure PDFs |
| `tools.lock.json`, `TOOLS.md` | The pinned tools and how each was built |
| `patches/` | The PHONI and Augmented PHONI fixes the benchmark binaries carry |
| `toy/` | A small end-to-end run of the whole workflow, on Rockfish |
| `data/` | Input download and preparation scripts; `data/locations/` maps logical input names to Rockfish paths (lab-internal) |
| `tests/` | Unit tests (`python -m pytest tests`) |

## The workflow

The freeze benchmark has three parts. A **manifest** (`ms2exp.manifest`) lists jobs and cells and
pins every input, tool binary included, by sha256 under a logical name. A **cell** is one tool
configuration on one read set, run three times cold, with run 1's whole output checked against
truth. The **runner** (`ms2exp.run.runner`) runs one job on either platform: it stages and verifies
the inputs, runs the cells in a seeded random order, checks them, and writes each cell's records and
a done marker as soon as the cell finishes, so an interrupted job resumes where it stopped. One
**launcher** per platform starts runners: `launchers/rockfish` for staging, dry runs and counting
runs, and `launchers/aws` for the timed runs. `ms2exp.collect` then gathers the cells into one
summary per cell group and checks them across tools.

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

Every number plotted is in `figures/draft/plotted-points.csv`, with its median query and load
times. The load-time table is in `figures/draft/CAPTIONS.md`. Fig rep is not regenerated from these
results (`results/README.md` says why).

## Toy run

`toy/toy_prep.sh TOY_DIR MS2_REPO MS2_COMMIT [--submodule-repo NAME URL ...]` builds the tools and
small mtb152 indexes at an ms2 commit, cuts small read sets, and writes a toy manifest whose cells
cover every path of the runner. `toy/toy_run.sh TOY_DIR` runs it and collects. Both need Rockfish,
since the toy's inputs are Rockfish paths. While ms2 is private, pass git bundles of the repository
and both submodules, made by `toy/make_source_bundles.sh`.

## Tools

`tools.lock.json` pins ms2 and its submodules, ropebwt3, MONI, PHONI and Augmented PHONI, and the
runtime libraries shipped with them; `TOOLS.md` says how each binary was built and checked.

## License

MIT (`LICENSE`). The patches in `patches/` are under their upstream projects' licenses.
