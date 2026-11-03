#!/usr/bin/env python3
"""Data for Fig Threads (thread-scaling appendix) from the threads-2885ee8-da1041e8 cell summaries.

Usage:
  threads_tikz.py [--results results/threads-2885ee8-da1041e8] [--out figures/tikz/data/fig-threads]

Writes <out>/<cohort>.dat with one row per configuration and thread count: speed in Mbases/s of
query time (median of three cold runs) and speedup over the same configuration at one thread.
The 60 Mb SMT check cell is left out. fig-threads.tex reads these files.
"""
import argparse
import glob
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
KEYS = {('ms2', 'ms'): 'lcp', ('ms2', 'psi'): 'psi', ('ms2', 'phi'): 'phi', ('ms2', 'dual'): 'dual',
        ('ropebwt3', None): 'rl1'}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--results', default=str(REPO/'results'/'threads-2885ee8-da1041e8'))
    ap.add_argument('--out', default=str(REPO/'figures'/'tikz'/'data'/'fig-threads'))
    args = ap.parse_args()
    results, out = Path(args.results), Path(args.out)
    rows = {}
    for path in glob.glob(str(results/'*'/'cells'/'*'/'summary.json')):
        s = json.loads(Path(path).read_text())
        if '-1gb-' not in s['read_set'] or s['status'] != 'checked':
            continue
        key = KEYS[(s['tool'], s.get('layout') if s['tool'] == 'ms2' else None)]
        cohort = s['read_set'].split('-ms-')[0]
        rows.setdefault(cohort, {}).setdefault(key, {})[s['threads']] = s['bases']/s['median_query_seconds']/1e6
    out.mkdir(parents=True, exist_ok=True)
    for cohort, by_key in rows.items():
        lines = ['key layout x y speedup']
        for key, by_t in sorted(by_key.items()):
            for t in sorted(by_t):
                lines.append(f'{key} {key} {t} {by_t[t]:.4f} {by_t[t]/by_t[1]:.3f}')
        (out/f'{cohort}.dat').write_text('\n'.join(lines) + '\n')
        print(cohort, len(lines) - 1, 'rows')


if __name__ == '__main__':
    main()
