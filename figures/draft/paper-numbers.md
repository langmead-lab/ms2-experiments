# Paper numbers

Written by `python -m ms2exp.figures.paper_numbers` from the committed results. Each number the paper quotes (papers/matching-stats at 567038f) is computed from the run records and rounded to the precision the paper prints. Status: agrees; edit (a known disagreement, with the paper edit below); FAIL (the paper and the data disagree, or a recorded edit no longer applies); unchecked (the data are outside this repository).


## Abstract

| Number | Paper | Computed | Status |
|---|---|---|---|
| ms2-min over PHONI, chr22 100 haplotypes (also Introduction, Results) | 35 to 48 | 35 to 48 | agrees |
| ms2-min over ropebwt3 -l1, human reads, at least (also Introduction) | 22 | min 22.93 | agrees |
| ms2-dual over ropebwt3 -l1, human reads, at least | 19 | min 19.04 | agrees |
| ms2-psi over ropebwt3 -l1, human reads, at least | 15 | min 15.26 | agrees |
| ropebwt3 index as a percentage of ms2-min's, about (also Introduction) | 15 | 14.46 | agrees |
| ms2-psi and ms2-dual over ropebwt3 -s8, N=100 (also Results) | 25 to 34 | 25 to 34 | agrees |
| PromethION flow cells per ms2-min thread (also Introduction, Results) | 21 | 21 | agrees |

## Introduction

| Number | Paper | Computed | Status |
|---|---|---|---|
| ms2-psi and ms2-dual over ropebwt3 -s6, N=100 (also Results) | 8.8 to 11 | 8.8 to 11 | agrees |

## Results: Recovery cost for ms2 modes

| Number | Paper | Computed | Status |
|---|---|---|---|
| mtb152 genome pieces, Mbases | 59.1 | 59.1 | agrees |
| HPRC2 genome pieces, Mbases | 59.9 | 59.9 | agrees |
| genome pieces sampled per index | 600 | 600 to 600 | agrees |
| simulated HiFi reads per index, Mbases | 60 | 60 to 60 | agrees |
| simulated HiFi median read length, kbases | 20 | 20 to 20 | agrees |
| recoveries per base, mtb152 genomes (also caption) | 0.00065 | 0.00065 | agrees |
| recoveries per base, HPRC2 genomes (also caption) | 0.00055 | 0.00055 | agrees |
| recoveries per base, mtb152 HiFi reads (also caption) | 0.021 | 0.021 | agrees |
| recoveries per base, HPRC2 HiFi reads (also caption) | 0.036 | 0.036 | agrees |
| ms2-min and ms2-dual above their oracles on HPRC2 genomes, within % | 6 | max 5.19 | agrees |
| HPRC2 genomes, ms2-min ns per base | 42 | 42 | agrees |
| HPRC2 genomes, ms2-dual ns per base | 40 | 40 | agrees |
| HPRC2 genomes, ms2-min oracle ns per base | 42 | 42 | agrees |
| HPRC2 genomes, ms2-dual oracle ns per base | 38 | 38 | agrees |
| HPRC2 genomes, ms2-dual operations per base | 1.01 | 1.01 | agrees |
| HPRC2 genomes, ms2-psi operations per base | 3.41 | 3.41 | agrees |
| HPRC2 genomes, ms2-phi operations per base | 1.07 | 1.07 | agrees |
| HPRC2 genomes, ms2-psi ns per base | 68 | 68 | agrees |
| HPRC2 genomes, ms2-phi ns per base | 47 | 47 | agrees |
| mtb152 genomes, ms2-min ns per base | 31 | 31 | agrees |
| mtb152 genomes, ms2-dual ns per base | 33 | 33 | agrees |
| mtb152 genomes, ms2-phi ns per base | 36 | 36 | agrees |
| mtb152 genomes, ms2-psi ns per base | 47 | 47 | agrees |
| HiFi reads over genomes, recoveries per base (times) | 33 to 65 | 33 to 65 | agrees |
| HiFi reads, ms2-min above its oracle, within % | 9 | max 8.63 | agrees |
| HiFi reads, ms2-dual slower than its oracle (%) | 32 to 37 | 32 to 37 | agrees |
| HiFi reads on HPRC2, ms2-phi time over its oracle (times) | 5.9 | 5.9 | agrees |

## Results: Comparison on chr22

| Number | Paper | Computed | Status |
|---|---|---|---|
| simulated reads per type, Mbases | 60 | 60 to 60 | agrees |
| ms2-min Mbases/s, 100 haplotypes | 22 to 29 | 22 to 29 | agrees |
| ms2-min over the fastest of MONI, PHONI and AugPHONI | 32 to 35 | 32 to 35 | agrees |
| ms2-min peak GiB | 1.0 to 1.4 | 1.0 to 1.4 | agrees |
| MONI, PHONI and AugPHONI peak GiB | 0.9 to 1.5 | 0.9 to 1.5 | agrees |
| ms2-phi Mbases/s | 4.7 to 12.6 | 4.7 to 12.6 | agrees |
| ms2-phi over the fastest of MONI, PHONI and AugPHONI | 7.7 to 14 | 7.7 to 14 | agrees |
| ropebwt3 peak GiB | 0.2 | 0.2 | agrees |
| ms2-min over ropebwt3 -l1 (times) | 14 to 32 | 14 to 32 | agrees |
| speed change of ms2-min, ms2-psi and ms2-dual from 100 to 10 haplotypes, within % | 6 | max 6.07 | agrees |
| ms2-phi on 10 haplotypes over 100 (times) | 1.6 to 2.3 | 1.6 to 2.3 | agrees |
| ms2-phi Mbases/s on 10 haplotypes | 11 to 20 | 11 to 20 | agrees |

## Results: Reads

| Number | Paper | Computed | Status |
|---|---|---|---|
| human_short reads | 480,000 | 480,000 | agrees |
| human_short Mbases | 60.0 | 60.0 | agrees |
| human_pacbio reads | 2,869 | 2,869 | agrees |
| human_pacbio Mbases | 60.0 | 60.0 | agrees |
| human_ont_modern reads | 2,706 | 2,706 | agrees |
| human_ont_modern Mbases | 59.9 | 59.9 | agrees |
| guccione_skin reads | 438,969 | 438,969 | agrees |
| guccione_skin Mbases | 60.0 | 60.0 | agrees |
| every 60 Mbase set, Mbases | 59.9 to 60.0 | 59.9 to 60.0 | agrees |
| PacBio median read length, kbases | 20.1 | 20.1 | agrees |
| Nanopore shortest read, bp | 404 | 404 | agrees |
| Nanopore longest read, kbases | 417 | 417 | agrees |
| Nanopore median read length, kbases | 3.4 | 3.4 | agrees |
| Nanopore N50, kbases | 91.7 | 91.7 | agrees |

## Results: Matching statistics on real reads

| Number | Paper | Computed | Status |
|---|---|---|---|
| oracle Mbases/s, every mode and set | 17 to 31 | 17 to 31 | agrees |
| human reads, ms2-min slower than its oracle (times) | 1.04 to 1.3 | 1.04 to 1.3 | agrees |
| human reads, ms2-psi and ms2-dual slower than their oracles (times) | 1.2 to 2.7 | 1.2 to 2.7 | agrees |
| skin, ms2-min Mbases/s | 10.2 | 10.2 | agrees |
| skin, ms2-min oracle Mbases/s | 26.1 | 26.1 | agrees |
| skin, ms2-psi Mbases/s | 2.8 | 2.8 | agrees |
| skin, ms2-dual Mbases/s | 5.2 | 5.2 | agrees |
| skin, ms2-psi and ms2-dual oracles Mbases/s, about | 18 | 18.11 | agrees |
| ms2-phi Mbases/s | 0.45 to 7.1 | 0.45 to 7.1 | agrees |
| ms2-phi as a percentage of ms2-psi's speed | 16 to 41 | 16 to 41 | agrees |
| ms2-min over ropebwt3 -l1 (times) | 23 to 130 | 23 to 130 | agrees |
| ropebwt3 -l1 Mbases/s | 0.08 to 1.12 | 0.08 to 1.12 | agrees |
| ms2-phi over ropebwt3 -l1 (times) | 2.6 to 6.3 | 2.6 to 6.3 | agrees |
| ms2-min Mbases/s on the Nanopore reads | 23.5 | 23.5 | agrees |

## Results: Locating SMEM occurrences

| Number | Paper | Computed | Status |
|---|---|---|---|
| human_short prefix reads | 48,000 | 48,000 | agrees |
| human_short prefix Mbases | 6 | 6 | agrees |
| human_pacbio prefix reads | 953 | 953 | agrees |
| human_pacbio prefix Mbases | 20 | 20 | agrees |
| human_ont_modern prefix reads | 494 | 494 | agrees |
| human_ont_modern prefix Mbases | 10 | 10 | agrees |
| guccione_skin prefix reads | 3,698 | 3,698 | agrees |
| guccione_skin prefix Mbases | 0.5 | 0.5 | agrees |
| ms2 locate index resident GiB | 337 | 337 | agrees |
| ms2-psi and ms2-dual million locations/s | 5.4 to 8.7 | 5.4 to 8.7 | agrees |
| ropebwt3 -s8 GiB | 97 | 97 | agrees |
| ropebwt3 -s6 GiB | 345 | 345 | agrees |
| ropebwt3 -s11 GiB | 25 | 25 | agrees |
| ms2-phi million locations/s | 4.0 to 8.5 | 4.0 to 8.5 | agrees |
| short reads, ms2-phi over ms2-dual | 0.60 | 0.60 | agrees |
| other sets, ms2-phi over ms2-dual | 0.79 to 0.98 | 0.79 to 0.98 | agrees |
| ropebwt3 -s11 thousand locations/s | 25 to 30 | 25 to 30 | agrees |

## Results: Thread scaling

| Number | Paper | Computed | Status |
|---|---|---|---|
| ms2 modes, 32 threads over 1 | 21.5 to 27.6 | 21.5 to 27.6 | agrees |
| ropebwt3 -l1, 32 threads over 1 | 25 to 30 | 25 to 30 | agrees |
| ms2-min Mbases/s at 32 threads | 469 to 498 | 469 to 498 | agrees |
| ms2-phi and ropebwt3, 64 threads over 32 | 1.16 to 1.28 | 1.16 to 1.28 | agrees |
| ms2-min and ms2-dual, 64 threads over 32, at most ("did not" benefit) | 1.00 | max 0.99 | agrees |

## Appendix: Fig threads caption

| Number | Paper | Computed | Status |
|---|---|---|---|
| short-read superset, reads | 7,982,098 | 7,982,098 | agrees |
| short-read superset, Mbases | 997.8 | 997.8 | agrees |
| Nanopore superset, reads | 43,346 | 43,346 | agrees |
| Nanopore superset, Mbases | 1,000.0 | 1,000.0 | agrees |

## Appendix: Fig rep caption

| Number | Paper | Computed | Status |
|---|---|---|---|
| (A) ms2-psi mean, exact reads | 1,471 | data outside the repository | unchecked |
| (A) ms2-psi mean at 0.14% errors | 41 | data outside the repository | unchecked |
| (A) ms2-psi mean, level near | 22 | data outside the repository | unchecked |
| (A) ms2-phi mean | 35 to 44 | data outside the repository | unchecked |
| (B) ms2-phi mean, 5 to 464 haplotypes | 2.4 to 143 | data outside the repository | unchecked |
| (B) ms2-psi mean | 35 to 45 | data outside the repository | unchecked |
| (B) ms2-dual mean | 0.7 to 2.5 | data outside the repository | unchecked |

## Appendix: Table A2 (load times)

| Number | Paper | Computed | Status |
|---|---|---|---|
| ms2-min, matching statistics: resident GiB | 100 | 100 | agrees |
| ms2-min, matching statistics: load s | 41.9 | 41.9 | agrees |
| ms2-min, matching statistics: GiB/s | 2.4 | 2.4 | agrees |
| ms2-psi, matching statistics: resident GiB | 97 | 97 | agrees |
| ms2-psi, matching statistics: load s | 45.2 | 45.2 | agrees |
| ms2-psi, matching statistics: GiB/s | 2.1 | 2.1 | agrees |
| ms2-phi, matching statistics: resident GiB | 147 | 147 | agrees |
| ms2-phi, matching statistics: load s | 60.7 | 60.7 | agrees |
| ms2-phi, matching statistics: GiB/s | 2.4 | 2.4 | agrees |
| ms2-dual, matching statistics: resident GiB | 228 | 228 | agrees |
| ms2-dual, matching statistics: load s | 94.7 | 94.7 | agrees |
| ms2-dual, matching statistics: GiB/s | 2.4 | 2.4 | agrees |
| ms2-psi, ms2-phi, ms2-dual, locating: resident GiB | 337 | 337 | agrees |
| ms2-psi, ms2-phi, ms2-dual, locating: load s | 186.7 | 186.7 | agrees |
| ms2-psi, ms2-phi, ms2-dual, locating: GiB/s | 1.8 | 1.8 | agrees |
| ropebwt3 -l1, matching statistics: resident GiB | 15 | 15 | agrees |
| ropebwt3 -l1, matching statistics: load s | 4.2 | 4.2 | agrees |
| ropebwt3 -l1, matching statistics: GiB/s | 3.4 | 3.4 | agrees |
| ropebwt3 -s11, locating: resident GiB | 25 | 25 | agrees |
| ropebwt3 -s11, locating: load s | 7.2 | 7.2 | agrees |
| ropebwt3 -s11, locating: GiB/s | 3.5 | 3.5 | agrees |
| ropebwt3 -s8, locating: resident GiB | 97 | 97 | agrees |
| ropebwt3 -s8, locating: load s | 27.6 | 27.6 | agrees |
| ropebwt3 -s8, locating: GiB/s | 3.5 | 3.5 | agrees |
| ropebwt3 -s6, locating: resident GiB | 345 | 345 | agrees |
| ropebwt3 -s6, locating: load s | 98.6 | 98.6 | agrees |
| ropebwt3 -s6, locating: GiB/s | 3.5 | 3.5 | agrees |

## Appendix: Index loading time

| Number | Paper | Computed | Status |
|---|---|---|---|
| ms2 load rates, GiB/s | 1.8 to 2.4 | 1.8 to 2.4 | agrees |
| ropebwt3 load rate, GiB/s, about | 3.5 | 3.48 | agrees |
| short reads, load over query time, ms2-min, ms2-psi and ms2-dual | 9 to 22 | 9 to 22 | agrees |
| short reads, load over query time, ms2-phi | 2 | 2 | agrees |

## Appendix: Measurement setup

| Number | Paper | Computed | Status |
|---|---|---|---|
| THP always over the default, ms2 on HPRC2 (times) | 1.4 to 2.2 | 1.4 to 2.2 | agrees |
| THP, ms2 matching statistics median | 1.8 | 1.8 | agrees |
| THP, ms2 locating median | 1.6 | 1.6 | agrees |
| THP, ropebwt3 median, about | 1.2 | 1.17 | agrees |
| THP, change in peak memory, at most % ("without changing peak memory") | 0 | max 0.03 | agrees |

## Appendix: Chromosome 22 with 10 haplotypes

| Number | Paper | Computed | Status |
|---|---|---|---|
| ms2-min Mbases/s on 10 haplotypes | 23 to 30 | 23 to 30 | agrees |
| ms2-min change from 100 haplotypes, within % | 6 | max 5.28 | agrees |
| MONI, PHONI and AugPHONI Mbases/s | 0.5 to 1.0 | 0.5 to 1.0 | agrees |
| ms2-min peak GiB | 0.6 to 1.0 | 0.6 to 1.0 | agrees |
| MONI, PHONI and AugPHONI peak GiB | 0.6 to 1.0 | 0.6 to 1.0 | agrees |

## Appendix: Locating with N=1 and N=10

| Number | Paper | Computed | Status |
|---|---|---|---|
| ms2-psi and ms2-dual over ropebwt3 -s6 at N=1 and N=10 | 5 to 11 | 5 to 11 | agrees |
| ms2-psi and ms2-dual thousand locations/s at N=1 | 111 to 364 | 111 to 364 | agrees |
| ms2-psi and ms2-dual million locations/s at N=10 | 0.92 to 2.52 | 0.92 to 2.52 | agrees |
| short reads, ms2-phi over ms2-dual at N=1 | 0.43 | 0.43 | agrees |
| short reads, ms2-phi over ms2-dual at N=10 | 0.45 | 0.45 | agrees |
| other sets, ms2-phi over ms2-dual at N=1 and N=10 | 0.68 to 0.90 | 0.68 to 0.90 | agrees |

## Not checked here

- Results: Reads: PacBio predicted mean error rate, "0.18%". Computed from the FASTQ base qualities, which the FASTA read sets do not keep.
- Results: Recovery cost for ms2 modes: simulated HiFi mean predicted accuracy, "99.8%". In the read sets' simulation manifests on Rockfish, not in data/read-sets.json.
- Results: Comparison on chr22: APR008 chr22 long arm, "34.0 Mbases". In the source haplotype record on Rockfish (d13_prep.py apr008-chr22), not in this repository.
- Results: Reads: Nanopore sampling rate and seed, "0.8%, seed 20261003". A preparation setting (data/download/prepare_hg002.py), not a measurement.
- Appendix: Data availability: HPRC2 precursor builds, "147 GiB, about 13 hours, 48 threads, 171 GiB; 83 and 331 GiB, about 10 hours, 46 to 48 threads". Build records of Rockfish runs outside the freeze manifests.

## Notes

- THP comparison: 279 cells in both runs, 230 timed three times in each (49 single-run count and Fig rep cells were copied into the rerun, identical); 142 of the timed cells are on HPRC2 (ms2 MS 56, ms2 locate 30, ropebwt3 MS 20, ropebwt3 locate 36); ropebwt3 speedups 1.09 to 1.26.

## Paper edits outstanding

None.

Totals: 148 agrees, 0 edit, 7 unchecked, 0 FAIL.
