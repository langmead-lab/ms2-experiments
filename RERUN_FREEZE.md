# Rerunning the freeze

The "freeze" is the recorded benchmark run behind every timed number in the paper. Its records are in `results/freeze-2885ee8-00e665de/`, and regenerating the figures from them takes a minute on a laptop (see "Figures" in `README.md`). This document is for a reader who wants to run the benchmarks again from public sources. It puts the existing pieces in order: the inputs, the job manifest, the AWS launcher, and the figures.

Before starting, run the small "toy" run in `README.md`. It goes through the same runner, checks and collection on one machine in about 40 minutes, and needs no AWS account.

The steps below have been checked piece by piece, not as one end-to-end run from public sources. The manifest step was checked to give the same 83 jobs and 294 cells as the paper's run. The input rebuild was checked on Rockfish (`evidence/hprc2-rebuild-20261107/`). The launcher is the one that produced the paper's records.

## Terms

- **Input**: any file a benchmark job reads (a read set, its truth, an index precursor, a tool binary), named by a logical path such as `precursors/human472/min` and pinned by sha256 in `data/inputs.tsv`.
- **Precursors**: the files ms2 builds its indexes from: the RLBWT runs (`bwt.heads`, `bwt.len`), LCP minima (`min`) and phi samples (`phi`).
- **Manifest**: the JSON list of jobs. A **job** is one index group and one read set; each job lists its **cells**, one tool configuration each. A manifest has an id such as `freeze-2885ee8-00e665de`.
- **Box**: one AWS EC2 instance running jobs one at a time.
- **SSA**: ropebwt3's sampled suffix array, used to locate matches; `-s8` and `-s6` are sampling intervals.

## What it takes

| Need | Amount |
|---|---|
| Machine for the input rebuild | Linux x86_64 with 48 threads, at least 171 GiB of memory for TeraLCP, and about 600 GiB of disk plus TeraLCP's temporary files |
| Time for the input rebuild | about 33 hours if run one after another: TeraLCP about 13 h, then two SSAs about 10 h each (the SSAs can run at the same time on two machines) |
| Downloads | about 46 GiB from Zenodo and GitHub releases (read sets, chr22 and mtb152 precursors, legacy indexes, tool binaries, Heng Li's human472 record) |
| Inputs uploaded to S3 | about 621 GiB (134 inputs), plus the indexes the build jobs write there (about 1.6 TiB) |
| AWS machines | i4i.metal spot instances in one region: 2 sockets of 32 cores, 1 TiB of memory, 8 local NVMe disks. The HPRC2 locate index alone is 337 GiB resident. The paper's runs used 5 boxes at a time |
| AWS account setup | an S3 bucket, an existing EC2 instance role, a spot vCPU quota of 128 per box (640 for 5 boxes), and the AWS CDK CLI |
| AWS instance-hours and cost of a full run | TBD. The design record's estimate before the runs was about 110 to 125 instance-hours (`docs/BENCH_FREEZE_DESIGN.md`, "First-pass AWS cost estimate"); the paper's runs were not costed afterwards |
| S3 storage | TBD; the design record estimated about $25 to $40 per month for the inputs and indexes |

## 1. Get the code and tools

```bash
git clone https://github.com/langmead-lab/ms2-experiments.git && cd ms2-experiments
export MS2EXP_DATA=/big/disk/ms2exp-data       # every input lands here under its logical name
python3 -m ms2exp.fetch tools/ropebwt3-cfad8475/ropebwt3 tools/runtime/libgcc_s.so.1 tools/runtime/libgomp.so.1
```

The fetch downloads the pinned ropebwt3 binary from the `tools-v1` release, which the input rebuild uses.

## 2. Rebuild the HPRC2 inputs that are not deposited

Seven large HPRC2 inputs (575 GiB) are not on Zenodo: the human472 FMD, its four precursors and its `-s8` and `-s6` SSAs. `data/prep/rebuild_hprc2.sh` makes them from Heng Li's public `human472.fmr` and checks each against its sha256. Its header lists the stages and their costs.

```bash
data/prep/rebuild_hprc2.sh /big/disk/hprc2-work tools fmd runs lcp     # about 13 h, 171 GiB of memory
data/prep/rebuild_hprc2.sh /big/disk/hprc2-work ssa                    # about 20 h for both SSAs
```

## 3. Fetch every other input

The paper's manifest names every input. This writes the names of those with a fixed sha256 and fetches them. Inputs from step 2 are verified and kept; the rest download from the three Zenodo deposits, Li's record and the `tools-v1` release.

```bash
python3 -c "import json; m = json.load(open('results/freeze-2885ee8-00e665de/manifest.json')); \
print('\n'.join(n for n, v in sorted(m['inputs'].items()) if len(v['sha256']) == 64 and not n.startswith('source/')))" > freeze-names.txt
python3 -m ms2exp.fetch --names freeze-names.txt
```

## 4. Write the manifest

`ms2exp.catalog from-results` takes the checksum of every input from the paper's run, and `ms2exp.manifest` writes a new manifest with the same jobs and cells. `--all-aws` adds the build jobs, which build ms2 and every ms2 index on AWS from the precursors. `--ms2-source repo` builds ms2 from the public repository at the paper's commit. `--allow-pending` lets the oracle traces be recorded on AWS during the run.

```bash
python3 -m ms2exp.catalog from-results --results results/freeze-2885ee8-00e665de --out inputs.json
python3 -m ms2exp.manifest --matrix freeze --all-aws --ms2-source repo --allow-pending --inputs inputs.json \
    --ms2-repo https://github.com/langmead-lab/ms2.git --ms2-commit 2885ee86270b4187fc45dcc69efd7de48a4aac03 \
    -o manifest.json
```

It prints the new manifest's id, for example `freeze-2885ee8-23b334e8: 83 jobs, 294 cells, 208 inputs (32 pending)`. The 32 pending inputs are the oracle traces. The id differs from the paper's because the manifest records when it was made.

## 5. Run the jobs on AWS

Follow `launchers/aws/README.md`: configure the stack, upload the code, manifest and inputs with `publish.sh`, deploy, watch the `job.done` markers, and destroy the stack when every job is done. Results go to `results/MANIFEST_ID/` under the S3 root.

## 6. Download the results and regenerate the figures

```bash
ID=freeze-2885ee8-23b334e8                      # the id step 4 printed
aws s3 sync s3://YOUR-BUCKET/freeze/results/$ID results/$ID \
    --exclude '*.canonical.gz' --exclude '*.oracle' --exclude '*.replog'
cp manifest.json results/$ID/manifest.json
python -m venv .venv && .venv/bin/pip install -e '.[figures,test]'
.venv/bin/python -m ms2exp.figures --results results --manifest $ID --out figures-rerun
```

The excluded files are the large locate outputs, oracle traces and recovery logs (about 24 GB in the paper's run); the figures do not need them. `figures-rerun/draft/paper-numbers.md` then compares every number the paper quotes with the new run.

Two appendix figures come partly from other runs, which the figures driver reads from `results/` unless told otherwise:

- Fig. A2 (operations per recovery) adds four batches counted on Rockfish after the freeze (`results/figrep-rungs-2885ee8/`).
- Fig. A3 (thread scaling) comes from its own manifest, `--matrix threads` (`docs/THREAD_SCALING_PLAN.md`, `results/threads-2885ee8-da1041e8/`), run with the same launcher.

## What will differ

Timings will not match the paper's to the digit: they depend on the box, the spot market's placement and the kernel. The checks that should match exactly are each cell's status, its truth check, and its operation and location counts. `ms2exp.collect` reports those for every cell group; `toy/compare.py` makes the same comparison, without timings, for the toy.
