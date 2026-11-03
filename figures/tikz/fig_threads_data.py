#!/usr/bin/env python3
"""Data for Fig Threads (thread-scaling appendix) from the threads-2885ee8-da1041e8 cell summaries.

Writes data/fig-threads/<cohort>.dat with one row per configuration and thread count: speed in
Mbases/s of query time (median of three cold runs) and speedup over the same configuration at
one thread. The 60 Mb SMT check cell is left out. Run from figures/tikz.
"""
import glob
import json
from pathlib import Path

RESULTS = Path(__file__).resolve().parents[2]/'results'/'threads-2885ee8-da1041e8'
OUT = Path(__file__).resolve().parent/'data'/'fig-threads'
KEYS = {('ms2', 'ms'): 'lcp', ('ms2', 'psi'): 'psi', ('ms2', 'phi'): 'phi', ('ms2', 'dual'): 'dual',
        ('ropebwt3', None): 'rl1'}


def main():
    rows = {}
    for path in glob.glob(str(RESULTS/'*'/'cells'/'*'/'summary.json')):
        s = json.loads(Path(path).read_text())
        if '-1gb-' not in s['read_set'] or s['status'] != 'checked':
            continue
        key = KEYS[(s['tool'], s.get('layout') if s['tool'] == 'ms2' else None)]
        cohort = s['read_set'].split('-ms-')[0]
        rows.setdefault(cohort, {}).setdefault(key, {})[s['threads']] = s['bases']/s['median_query_seconds']/1e6
    OUT.mkdir(parents=True, exist_ok=True)
    for cohort, by_key in rows.items():
        lines = ['key layout x y speedup']
        for key, by_t in sorted(by_key.items()):
            for t in sorted(by_t):
                lines.append(f'{key} {key} {t} {by_t[t]:.4f} {by_t[t]/by_t[1]:.3f}')
        (OUT/f'{cohort}.dat').write_text('\n'.join(lines) + '\n')
        print(cohort, len(lines) - 1, 'rows')


if __name__ == '__main__':
    main()
