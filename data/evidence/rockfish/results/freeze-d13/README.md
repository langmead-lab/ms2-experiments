# Freeze inputs for design decision D13 (Oct 30 2026)

Prepared by `bench/freeze/d13_prep.py`; every set directory under `inputs/` has a `manifest.json`
with its source, rule, counts and checksums. Logical names and hashes are in
`freeze-catalog/inputs-d13.json`.

- chr22 simulated sets (Fig 4, A2, A1 series B): simulated from the long arm of chromosome 22 of
  APR008 haplotype 1 (GCA_050494745.1, scaffold JBHMUJ010000002.1, 34.0 Mb, no N bases).
  minimap2 2.28 -x asm5 of all 79 scaffolds against CHM13 v2.0 chr22 (job 32742987,
  `apr008-chr22/`) shows this scaffold covers chr22 17.4 to 51.3 Mb (33.85 Mb in alignments of at
  least 500 kb, 99.5 percent of the scaffold) and no other scaffold aligns there. The short arm
  and centromere (0 to 17.4 Mb) are split over at least 10 scaffolds shared with the other
  acrocentrics and were excluded. Approved by Ben Langmead ("Yes, long arm"). Simulators,
  models, settings and seeds are chr22_sim.py's; truth is ropebwt3 cfad8475 mem -t1 -l1 on the
  chr22-h10 and chr22-h100 FMDs (job 32746173).
- Every 60 Mb set: ACGT-only whole reads first, then whole reads in sampled order up to at most
  60,000,000 bases. Short and skin were rebuilt from larger pools (short: the real-hprc2 core,
  truth sliced; skin: ERR14066769, new truth on the human472 FMD, job 32743099).
- Fig 3 pieces: seeded samples of 600 records (seeds 2026103001 HPRC2, 2026103002 mtb152).
- Locate sets: whole-read prefixes of the 60 Mb sets (6, 20, 10 and 0.5 Mb).
