# Expected toy results

The collected group summaries, collect report and manifest of the toy run made on 2026-11-03 on
devlangmead1 (the lab's development server on Rockfish) with `toy_prep.sh` and `toy_run.sh` from
this repository, ms2 2885ee86 (`paper-v1`) from GitHub, and inputs linked from the lab's copies by
`python -m ms2exp.fetch --set toy --link`. All 34 cells finished: 33 checked against truth, and
the `phi-ordinary-capped` cell capped by design (its cap is 1 s), and every locate comparison with
ropebwt3 agreed. `toy/compare.py --expected toy/expected TOY_DIR` checks another run against these
files. Timing fields are devlangmead1's, from runs without root (inputs evicted with
posix_fadvise), so they are for orientation only.
