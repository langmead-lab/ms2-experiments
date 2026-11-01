#!/usr/bin/env python3
"""Violin densities of LCE steps per reposition for Fig Rep, from the ms-replog logs.

For each run and mode (psi, phi, dual), takes every reposition with a nonzero match length
(as fig_rep.py summarize does), and estimates the density of log10(steps) with a Gaussian
kernel of KERNEL_DECADES on a fine histogram. The density is sampled at POINTS heights
between the 0.5th and 99.5th percentiles and scaled so that its peak is 1. Besides the three modes,
'min' is the per-repositioning minimum of the psi and phi costs. Writes one JSON object:
{run: {mode: {"log10_steps": [...], "density": [...], "mean": m, "n": n}}}.

Usage (on Rockfish, in results/campaign-20261012/fig_rep): fig_rep_violin.py [--plus-one] OUT.json RUN_DIR...
"""

import json
import sys
from pathlib import Path

import numpy as np

from ms2exp.figures.fig_rep import MODES, numpy_dtype

KERNEL_DECADES = 0.06
BINS_PER_DECADE = 200
POINTS = 90


def density(steps):
    x = np.log10(steps)
    lo, hi = np.percentile(x, [0.5, 99.5])
    edges = np.arange(np.floor(x.min()) - 0.5, np.ceil(x.max()) + 0.5, 1 / BINS_PER_DECADE)
    hist, _ = np.histogram(x, bins=edges)
    sigma = KERNEL_DECADES * BINS_PER_DECADE
    k = np.arange(-int(4 * sigma), int(4 * sigma) + 1)
    smooth = np.convolve(hist, np.exp(-0.5 * (k / sigma) ** 2), mode='same')
    centers = (edges[:-1] + edges[1:]) / 2
    # A degenerate spread (all values equal) still gets a short violin.
    if hi - lo < 0.02:
        lo, hi = lo - 0.01, hi + 0.01
    grid = np.linspace(lo, hi, POINTS)
    d = np.interp(grid, centers, smooth)
    return grid, d / d.max()


def main():
    out, dt = {}, numpy_dtype()
    # --plus-one shifts every cost by one step, so that repositions settled without a walk (zero
    # steps, common in ms2-phi and ms2-dual) stay on the log axis.
    args = sys.argv[1:]
    plus_one = '--plus-one' in args
    args = [a for a in args if a != '--plus-one']
    shift = 1.0 if plus_one else 0.0
    for run in args[1:]:
        run = Path(run)
        logs = {m: np.fromfile(run/f'{m}.replog', dtype=dt) for m in MODES}
        walked = logs['psi']['cap'] > 0
        out[run.name] = {}
        costs = {m: (logs[m]['psi'] + logs[m]['phi'])[walked].astype(np.float64) + shift for m in MODES}
        # 'min' is not a method: for each repositioning, the smaller of ms2-psi's and ms2-phi's
        # own step counts, the cost of an oracle that knew in advance which walk would finish first.
        costs['min'] = np.minimum(costs['psi'], costs['phi'])
        for m in (*MODES, 'min'):
            steps = costs[m]
            grid, d = density(steps)
            out[run.name][m] = dict(log10_steps=[round(float(v), 4) for v in grid],
                                    density=[round(float(v), 4) for v in d],
                                    mean=float(steps.mean()), n=int(len(steps)))
        print(run.name, {m: round(out[run.name][m]['mean'], 1) for m in out[run.name]}, flush=True)
    Path(args[0]).write_text(json.dumps(dict(kernel_decades=KERNEL_DECADES, points=POINTS, plus_one=plus_one,
                                                 range='0.5th to 99.5th percentile', runs=out)) + '\n')


if __name__ == '__main__':
    main()
