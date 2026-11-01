# Expected toy results

The collected summaries and manifest of the toy run made on 2026-11-01 on devlangmead1 with the
freeze code as it ran on AWS (`freeze-code.tar.gz`, sha256 300e2338...), ms2 2885ee86 from source
bundles. A run of `toy_prep.sh` and `toy_run.sh` from this repository at the same commit gave the
same summaries, report and done markers apart from timings, memory, hosts and paths
(`compare.py`); the only other differences were the sha256s of the separately built mtb152
indexes, whose ORBITIDX headers record the build path and time. Timing fields here are
devlangmead1's, from runs without root (inputs evicted with posix_fadvise), so they are for
orientation only.
