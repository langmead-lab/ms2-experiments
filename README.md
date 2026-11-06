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
| `data/` | Where every input came from (`data/README.md`, the entry point), with the download and preparation scripts, input lists and evidence; `data/inputs.tsv` pins each input by sha256, and `data/locations/` maps logical input names to Rockfish paths (the Langmead lab's own copies) |
| `evidence/` | Records of the binary checks and byte-identical rebuilds that `TOOLS.md` cites |
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
| `ms2exp.read_sets` | The read-set records in `data/read-sets.json` |
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
times. The load-time table (Table A2) is `figures/tex/tab-load.tex`. `figures/draft/CAPTIONS.md`
holds draft captions built from the run records and `data/read-sets.json`. Every number the paper
quotes is listed, by id, in `figures/paper-values.json`; `python -m ms2exp.figures.paper_numbers`
computes each one from the records and writes `figures/draft/paper-numbers.md`, which also lists the
paper edits still outstanding (`--strict` fails on any disagreement). Fig rep is not regenerated
from these results (`results/README.md` says why).

## Small "toy" run

The toy runs every path of the runner (ms2 MS, oracle, counting and locate cells, ropebwt3, MONI,
PHONI and Augmented PHONI) on mtb152 and chr22-h10 inputs, then checks the results against
`toy/expected/`. It needs x86_64 Linux with glibc 2.28 or later, g++ with C++17 and static
libstdc++, make, zlib headers, git and Python 3.10 or later; numactl is optional. Expect about
5.8 GB of downloads (about 50 minutes at the 2 MB/s we saw from Zenodo), 15 GB of disk and under
1.5 GB of memory. On the lab's development server it
took about 11 minutes to prepare and 10 minutes to run, at 1.3 GB peak memory.

```bash
git clone https://github.com/langmead-lab/ms2-experiments.git && cd ms2-experiments
export MS2EXP_DATA=$HOME/ms2exp-data TOY=$HOME/ms2exp-toy
python3 -m ms2exp.fetch --set toy          # 34 inputs, each checked by sha256
toy/toy_prep.sh $TOY                       # builds ms2 at paper-v1, the toy indexes and manifest
toy/toy_run.sh $TOY                        # runs the 34 cells on this machine, then collects
python3 toy/compare.py --expected toy/expected $TOY
```

The fetch step downloads the read sets, precursors and legacy indexes from the three Zenodo
deposits ([10.5281/zenodo.23964218](https://doi.org/10.5281/zenodo.23964218),
[10.5281/zenodo.23964220](https://doi.org/10.5281/zenodo.23964220) and
[10.5281/zenodo.23964222](https://doi.org/10.5281/zenodo.23964222)); `data/inputs.tsv` gives each
input's sha256 and URL. The pinned ropebwt3, MONI, PHONI and Augmented PHONI binaries and the
three runtime libraries beside them (26 MB in all; `TOOLS.md`) are assets of this repository's
[`tools-v1` release](https://github.com/langmead-lab/ms2-experiments/releases/tag/tools-v1).
`python3 -m ms2exp.fetch --set toy --dry-run` lists every input with its size and source before
downloading anything. Inside the Langmead lab, `python3 -m ms2exp.fetch --set toy --link --from-locations
data/locations/*.json` links the lab's copies of every input instead.

`toy_run.sh` runs each job directly on this machine; it does not use Slurm. The protocol times runs
with transparent huge pages set to `always` (`sudo sh -c 'echo always >
/sys/kernel/mm/transparent_hugepage/enabled'`); without that setting the toy warns and runs anyway.
The last command prints `equivalent` when every status, check, operation count, location count,
command and input matches; timings are not compared. The toy builds ms2 with your compiler, so it
checks the workflow, not the bytes of the AWS freeze binaries; ropebwt3 and the legacy tools are the
pinned binaries. `toy/expected-v1-snapshot/` holds the earlier toy, on v1 read sets, that showed
`ms2exp` gives the same results as the code snapshot the freeze ran (`results/README.md`).

## Tools

`tools.lock.json` pins ms2 and its submodules, ropebwt3, MONI, PHONI and Augmented PHONI, and the
runtime libraries shipped with them; `TOOLS.md` says how each binary was built and checked.

## Acknowledgements

`ms2` builds on Orbit, the move-structure library by Nathaniel Brown, and TeraTools, the move-structure building library and its TeraMS component for computing matching statistics. Its psi and phi query engine follows the design and code of TeraTools (TeraMS and TeraLCP) by Ahsan Sanaullah, Shaojie Zhang and the S. Zhang Lab at UCF, and it builds from TeraLCP's output. fmd2rlbwt uses ropebwt3's FMD reader by Heng Li.

## License

Software in this repo is under the MIT license (`LICENSE`). The patches in `patches/` are under their upstream
projects' licenses.
