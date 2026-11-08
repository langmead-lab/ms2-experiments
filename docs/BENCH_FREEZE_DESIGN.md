# Freeze benchmark workflow design

Design record of the freeze benchmark, written October 29 and 30, 2026, before the runs. Decisions that changed during execution carry an Outcome note; the code and `results/README.md` describe what ran. Code docstrings and the committed result records cite this file by its bare name, `BENCH_FREEZE_DESIGN.md`, and its decisions as D1 to D13.

Goal: one workflow that regenerates the data for manuscript Figs 3, 4, 5, 6, A1 to A6 (and Table tab:load, which derives from the same runs), targetable at Rockfish (Slurm) or AWS (CDK, one metal instance type).

Figure map as planned: 3 ops, 4 ms22_100, 5 ms_oracle, 6 locate100, A1 rep, A2 ms22_10, A3 locate1, A4 locate10, A5 ms_ftab, A6 locate_ftab. Figure pipeline: paper_figures.py, paper_figures_tikz.py and fig_rep*.py, now `ms2exp.figures`.
Outcome: the decisions below keep these planned appendix numbers, which code comments and records also use. The paper numbers its figures differently: Figs. 3 to 6 as above, then A1 ms22_10, A2 rep, A3 threads (`THREAD_SCALING_PLAN.md`), A4 locate1 and A5 locate10. The two ftab figures (planned A5 and A6) are not in the paper. Table tab:load is the paper's Table A3, and Table tab:haplotypes its Table A2. The figure table in `README.md` maps every name.

## Decisions

### D1. Inputs and index reuse (Oct 29)
Reuse every artifact not produced by an ms2 command (read sets amended by D13): read sets and truth, source FASTAs, grlBWT output, TeraLCP minima, ropebwt3 FMD and SSAs (s6, s8, s11), MONI/PHONI/AugPHONI chr22 indexes. Re-run every build step that invokes ms2 with the final ms2 binary: build-rlbwt (ms2-min), tms-build psi/phi/full, tms-add-summary, index-add-seq-ids, ftab k=14, for HPRC2, chr22 h10/h100, mtb152, and the A1 chr22 ladder (D13: h1, h2, h5, h20, h50 beyond h10/h100).

### D2. Platform roles (Oct 29, provisional until cost estimates)
Rockfish is used for staging: index builds, a dry run of the full workflow, and correctness checks against truth. AWS metal produces the timed numbers that go in the paper.
Outcome: every job of the freeze, builds included, ran on AWS (see D3, D7 and D12).

Data volume facts: the HPRC2 ms2 build precursors are bwt.heads 5 GiB, bwt.len 25 GiB, min 30 GiB, records.tsv under 1 GiB (lcp_index 144 GiB is TeraLCP intermediate). Timed-run inputs total about 1.6 TiB, of which ms2 HPRC2 indexes are about 1.13 TiB (about 0.79 TiB if the seqids index serves MS dual too) and ropebwt3 FMD plus SSAs are 0.43 TiB. S3 Standard is about $0.023 per GB-month, so 1.75 TB is about $40/month; ingress and same-region S3 to EC2 transfer are free.

### Note: ms2 source location
The workflow takes the ms2 source as a parameter (repository URL plus commit) and records that commit in every result.

### Upload vs build time estimates (Oct 29)
Past Rockfish to S3 uploads with default aws-cli settings ran at about 88 to 100 MiB/s; tuned parallel uploads are untested.
- All timed-run inputs, 1.6 TiB: about 5 h at 95 MiB/s.
- ms2 HPRC2 indexes only, 0.79 to 1.13 TiB: about 2.5 to 3.5 h.
- ropebwt3 FMD plus SSAs, 0.43 TiB (needed in either option): about 1.3 h.
- ms2 build precursors, 60 GiB: about 11 min.
HPRC2 ms2 builds at e637ac0 on Rockfish: build-rlbwt 0.75 h, psi tms-build 0.5 h, full 10.3 h (369 GiB maxrss), phi 10.7 h (301 GiB maxrss); tms-add-summary, index-add-seq-ids and ftab times for the final binary are unknown. Critical path is about 11 to 13 h whether the builds run in parallel on Rockfish nodes or concurrently on one 1 TiB AWS box.

### D3. Where ms2 indexes are built (Oct 29)
Build every ms2 index once on Rockfish with the final ms2 commit, record sha256 checksums, run the dry run and truth checks against those files, then upload the same bytes to S3. AWS never builds indexes; the EC2 side only runs timed experiments.
Outcome: the freeze manifest is an all-AWS one (`ms2exp.manifest --all-aws`). Its build jobs built the ms2 binaries and every ms2 index on i4i.metal from the pinned precursors (`results/freeze-2885ee8-3c8c8967/build-*`), and the timed jobs read those outputs by sha256.

### D4. Instance type (Oct 29)
r6idn.metal in us-east-1 (Xeon Platinum 8375C Ice Lake, 2 sockets x 32 cores, 1024 GiB, 512 GiB per socket, 7.6 TB local NVMe), spot. Fallback pool: i4i.metal (same CPU, memory and topology). Oct 6 survey medians: r6idn.metal $1.42/h us-east-1; i4i.metal $1.71/h. Refresh before costing.
Outcome: the launcher listed r6idn.metal first, but spot capacity gave i4i.metal to every box; every committed record is i4i.metal.

### D5. One job per box (Oct 29)
Each EC2 box runs one experiment at a time, bound with numactl to socket 0 (cores and memory); socket 1 stays idle. Wall time is shortened, if needed, by running several boxes in parallel, which does not change cost.

### D6. Repetition protocol (Oct 29)
Three timed runs per cell, no warmup, median reported (matches the paper's "median of three timed runs"). Before every run, drop the page cache (echo 3 > /proc/sys/vm/drop_caches), so each index load is a cold read from local NVMe. The first run's output is checked against truth. Round order is shuffled across cells with a recorded seed. Table tab:load comes from these same runs. This replaces the exploratory scripts' one warmup plus two timed runs.

### Fig 3 inputs carried forward from FIG_OPS_LATENCY_EXPERIMENT.md
Genome query sets cut into 100 kb pieces (human472-apr008-hap1-p100k-v1, 30,424 pieces; mtb152-heldout16-p100k-v1, 718 pieces), --block-bytes 67108864 on all four Fig 3 sets, K = 32. The K = 1/4 and unchunked controls were diagnostic and are not part of the official run.

### D7. Counting runs on Rockfish (Oct 29)
Deterministic operation counts run on Rockfish during staging, not on AWS: Fig 3 right panel (MS_STATS/TMS_STATS build) and all of A1 (TMS_REP_LOG build; mtb152 error series and the chr22 ladder of D13). Counting builds use the same final ms2 commit as the timing build, differing only in compile flags. A1's ladder-only indexes (h1, h2, h5, h20, h50) never go to S3.
Outcome: the counting and A1 jobs ran on AWS too, as untimed jobs of the all-AWS manifest.

### D8. MS cell lists (Oct 29)
- Fig 3 (AWS timing): {mtb152, HPRC2} x {100 kb genome pieces, PacBio} x {min, psi, phi, dual} x {ordinary, oracle}.
- Fig 4 / A2: chr22 {h100, h10} x {sim short, PacBio, ONT}; ms2 min/psi/phi/dual ordinary; ropebwt3 -l1/19/25/30/40 with -K 1m; MONI, PHONI, AugPHONI with one-read baseline runs.
- Fig 5: HPRC2 x 4 real 60 Mb sets; ms2 min/psi/phi/dual x {ordinary, oracle}; ropebwt3 -l1/19/25/30/40. ms2-phi is included unconditionally.
- A5: HPRC2 short reads; ms2 min/psi/phi/dual x {ordinary, ftab k=14}. Needs a phi ftab.
- Dropped from the official run: oracle-scan, chr22 real-read cells, locate oracle replay.

### D9. Locate cell lists and indexes (Oct 29)
- Fig 6 / A3 / A4: HPRC2 x 4 locate read sets (short 49,152 reads, PacBio 20 Mb, ONT 10 Mb, skin 0.4 Mb) x caps {1, 10, 100}; ms2-psi on a new HPRC2 psi-locate index (branch psi-locate, ec43f8d: full layout without PHI_LEN); ms2-phi and ms2-dual on the shared full + seqids index; ropebwt3 -l1 -pN -K 1m at -s11, -s8, -s6.
- A6: short reads, caps {1, 10, 100}, ms2 psi/phi/dual x {ordinary, ftab}.
- The full + seqids file also serves every HPRC2 MS dual cell (Figs 3, 5, A5), replacing the plain tms.full file, provided the Rockfish dry run shows MS dual timing on it matches the plain file within noise.
- ms2-min does not locate.

## First-pass AWS cost estimate (Oct 29, from Rockfish e637ac0 run times; refine after the dry run)
| Group | Instance-hours |
|---|---|
| Fig 5 + A5 (adds phi) | about 12 |
| Fig 3 timing (mtb152 + HPRC2) | about 9 |
| Fig 4 + A2 (chr22) | about 7 |
| Fig 6, A3, A4, A6 (drops locate oracle, adds phi) | about 70 to 80 |
| Boot, S3 pull, retries | about 10 to 15 |
| Total | about 110 to 125 |
At r6idn.metal spot (about $1.42/h, Oct 6) this is about $160 to $180; on-demand (about $10.4/h, unverified) about $1,150 to $1,300. S3 storage about $25 to $40 per month.

### D10. Spot with per-cell checkpoints (Oct 29)
Spot instances. A job is a list of cells; after each cell's three runs, its outputs, summary and a done marker go to S3 at once. A job skips cells that already have a done marker, so an interruption loses only the cell in progress. The orchestrator relaunches interrupted jobs, falling back to i4i.metal if r6idn.metal capacity is short. The same resume logic serves Rockfish dry runs that hit a Slurm time limit.

### D11. Architecture (Oct 29)
Code lived in the lab's working tree and is now the `ms2exp` package of this repository.
1. Manifest: a committed, declarative list of jobs (one per index group and read set), each listing its cells; every input is a (logical name, sha256) pair; ms2 and ropebwt3 commits and the D6 protocol are pinned. The D8/D9 matrices are generated into it once.
2. Cell runner: run_job.py <job-id>, identical on both platforms. Resolves inputs from a local root or S3, stages to local NVMe with checksum verification, runs cells under the D5/D6 protocol, checks against truth, writes per-cell results and done markers (D10). Reuses measurement.py, ms_check.cpp and the offset comparators of the exploratory runs.
3. Launchers: Rockfish sbatch wrapper (inputs from /vast; no root, so no page-cache drop; dry-run timings are for planning only). AWS CDK stack adapted from an earlier lab profiling stack: r6idn.metal spot instead of c5.metal, no perf/likwid; user-data builds NVMe RAID0, installs pinned binaries from S3, pulls job IDs from an S3 queue, runs run_job.py, and shuts the box down when the queue is empty.
paper_figures.py reads the per-cell summaries with path changes only.

### D12. Full-size Rockfish dry run (Oct 29)
Every cell, full read sets, on parallel bigmem nodes (about 110 to 125 node-hours, 1 to 2 days wall). It gates the AWS run: every cell passes truth, and the seqids-vs-plain dual MS timing check (D9) passes. Oracle traces recorded here ship to S3 with the indexes (valid because index bytes are identical, D3). Per-cell times replace the first-pass AWS estimate and drive job packing. Results also serve as a fallback dataset.
Outcome: no full Rockfish dry run took place. The all-AWS manifest checked run 1 of every cell against truth on AWS itself, and the toy run (`toy/`) exercised every path of the runner on Rockfish beforehand. The D9 check ran on AWS and passed (`results/README.md`).

### D13. Dataset consolidation and round sizes (Oct 30)
This amends D1 (read sets are no longer all reused as-is) and D8/D9 (query sets).
1. One held-out human haplotype: APR008 hap1 (GCA_050494745.1) is the source for Fig 3 HPRC2 queries and for the chr22 simulated reads of Fig 4, A2 and A1 series B, replacing NA19240 paternal chr22. Re-simulate chr22 short (ART HS25, 125 bp), PacBio HiFi (pbsim3 + ccs, real PacBio length distribution) and ONT (pbsim3 ERRHMM-ONT-HQ, 96%) from APR008 chr22, 60 Mb each, and recompute truth with ropebwt3 -l1 on each chr22 index. Resolved Oct 30: only the long arm is assignable. Source is APR008 hap1 scaffold JBHMUJ010000002.1 (34.0 Mb, covering CHM13 chr22 17.4 to 51.3 Mb); the short arm and centromere are split across scaffolds shared with other acrocentrics and are excluded.
2. Fig 3 genome queries: a seeded random sample of 600 of the 100 kb pieces (60 Mb) for HPRC2 (from 30,424 pieces of APR008 hap1) and for mtb152 (from 718 pieces of the held-out genomes plus reverse complements). Truth is sliced from the existing whole-genome truth. Every MS query set is then 60 Mbases.
3. A1 series B ladder: chr22 h1, h2, h5, h10, h20, h50, h100 (first N haplotypes in Table tab:haplotypes order), replacing h1/2/3/5/8/13/21/34/55. h10 and h100 are the Fig 4/A2 indexes; five ladder-only builds remain.
4. Every 60 Mb read set is exactly 60.0 Mb of whole reads: apply the A/C/G/T-only filter first, then take whole reads in sampled order, stopping at the last read that keeps the total at or below 60,000,000 bases (so every set reports as 60.0 Mb). Applies to the four real HPRC2 sets and all simulated sets.
5. Locate read sets are prefixes of the 60 Mb sets: the first 6 Mb (short), 20 Mb (PacBio), 10 Mb (ONT) and 0.5 Mb (skin). Equal sizes may replace these if dry-run per-cell times allow.
Optional trim, not adopted yet: ropebwt3 MS modes -l1, -l19, -l40 only.

## Defaults
- Fig 3 time basis is query time including query read and output write, excluding index load, as the paper text states.
- MONI, PHONI, AugPHONI query time stays whole-process wall minus a one-read run; under D6 both runs are cold.
- Truth sets are reused (D1); new truth is needed only where cells change output, which none do.
- HPRC2 phi ordinary cells keep a 3600 s cap per run; a cell that hits it is reported as capped, not dropped.
- ropebwt3 commit stays at cfad8475 (the e637ac0-era build).
- Locate work is split into AWS jobs by read set and tool so no job exceeds about 12 h.

## Timeline after the final ms2 commit is tagged
1. Rockfish ms2 builds: about 1 day (phi and full in parallel, then tms-add-summary, seqids, ftabs, psi-locate; chr22, mtb152 and A1 ladder alongside).
2. Rockfish full dry run plus counting runs (Fig 3 ops panel, A1): 1 to 2 days.
3. Upload to S3: about 3 h (ms2 HPRC2 indexes and traces about 0.8 to 1.1 TiB, ropebwt3 0.43 TiB, small inputs).
4. AWS timed runs: about 1 day wall on 6 to 8 r6idn.metal spot boxes.
5. paper_figures.py and paper_figures_tikz.py, copy data into the paper repo.
