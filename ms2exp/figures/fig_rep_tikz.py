#!/usr/bin/env python3
"""Fig Rep TikZ data: box statistics of LCE steps per reposition, from fig_rep.py summaries.

Writes reports/campaign-20261012/figures-tikz/data/fig-rep/:
  a-boxes.tex, b-boxes.tex   one \\figRepBox{mode}{x}{p5}{p25}{p50}{p75}{p95}{mean} per mode
                             and batch, and one \\figRepRatio{x}{...} per batch
  a-ticks.tex, b-ticks.tex   \\def\\FigRepTicksA{...} and \\def\\FigRepTickLabelsA{...}
  batches.tsv                every plotted number, one row per series, batch and mode
  a-violins.tex, b-violins.tex  with --violins: one \\figRepViolin{mode}{polygon} per mode and
                             batch, and one \\figRepMeanLine{mode}{points} per mode
fig-rep.tex draws them. Series A batches are labeled by measured repositions per base and
the error rate; series B by haplotype count.

Usage: fig_rep_tikz.py --a-dir DIR --b-dir DIR [--out DIR] [--violins JSON]
DIR holds runs/<run>/summary.json and inputs/<run>/manifest.json as copied from Rockfish.
"""

import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1]/'reports/campaign-20261012/figures-tikz/data/fig-rep'
MODES = ('psi', 'phi', 'dual')
# Series A: (run, errors per base as a label).
SERIES_A = [('exact', '0'), ('e0.0001', '0.01\\%'), ('e0.0002', '0.02\\%'), ('e0.0005', '0.05\\%'),
            ('e0.001', '0.1\\%'), ('e0.003', '0.3\\%'), ('e0.01', '1\\%'), ('e0.03', '3\\%'),
            ('e0.05', '5\\%'), ('e0.10', '10\\%'), ('e0.15', '15\\%')]
SERIES_B = ['h1', 'h2', 'h5', 'h10', 'h20', 'h50', 'h100']
# The redesigned sweeps (October 28) step each series' relevant mean by a near-constant
# factor: series A by error rates interpolated from the first sweep, series B by
# near-Fibonacci haplotype counts. The 0.5 percent batch extends A past psi's floor.
SERIES_A2 = [('exact', '0'), ('e0.000032', '0.003\\%'), ('e0.00011', '0.011\\%'), ('e0.00019', '0.019\\%'),
             ('e0.00030', '0.03\\%'), ('e0.00047', '0.047\\%'), ('e0.00076', '0.076\\%'), ('e0.0014', '0.14\\%'),
             ('e0.0030', '0.3\\%'), ('e0.0050', '0.5\\%')]
SERIES_B2 = ['h1', 'h2', 'h3', 'h5', 'h8', 'h13', 'h21', 'h34', 'h55']
# Ben's suggested haplotype counts, drawn for comparison.
SERIES_B3 = ['h1', 'h2', 'h3', 'h4', 'h5', 'h10', 'h25', 'h50']
# The freeze x axes redrawn (November 1): series A at every other freeze error rate from 0.011
# percent, plus 1 and 2 percent, which show psi's floor; series B at doubling haplotype counts
# from 5, plus 200 and every r2_assemblies_only haplotype (464) built for this design.
SERIES_A5 = [('exact', '0'), ('e0.00011', '0.011\\%'), ('e0.00030', '0.03\\%'), ('e0.00076', '0.076\\%'),
             ('e0.0014', '0.14\\%'), ('e0.0030', '0.3\\%'), ('e0.0050', '0.5\\%'), ('e0.010', '1\\%'),
             ('e0.020', '2\\%')]
SERIES_B5 = ['h5', 'h10', 'h20', 'h50', 'h100', 'h200', 'h464']
# Design: (A run prefix, A batches, B cells, file tags, macro suffixes).
DESIGNS = {'1': ('mtb152-A-', SERIES_A, SERIES_B, ('a', 'b'), ('A', 'B')),
           '2': ('mtb152-A2-', SERIES_A2, SERIES_B2, ('a2', 'b2'), ('Atwo', 'Btwo')),
           '3': ('mtb152-A2-', SERIES_A2, SERIES_B3, ('a2', 'b3'), ('Atwo', 'Bthree')),
           # The freeze runs: series A as redesigned, series B at the haplotype counts of the chr22
           # figures, written under the second design's names so that fig-rep.tex reads them unchanged.
           '4': ('mtb152-A2-', SERIES_A2, SERIES_B, ('a2', 'b2'), ('Atwo', 'Btwo')),
           # Also under the second design's names, for the same reason.
           '5': ('mtb152-A2-', SERIES_A5, SERIES_B5, ('a2', 'b2'), ('Atwo', 'Btwo'))}


def fmt(v):
    return f'{v:.6g}'


def rate(r):
    """Repositions per base to two significant digits, as plain text."""
    return f'{r:.2g}' if r >= 0.001 else f'{r:.1g}'


def write(out, tag, rows, tsv, T=None):
    lines = [f'% Fig Rep series {tag.upper()}: written by bench/campaign-20261012/fig_rep_tikz.py.']
    ticks, labels = [], []
    for x, (run, label, s) in enumerate(rows, start=1):
        ticks.append(str(x))
        labels.append('{' + label + '}')
        for mode in MODES:
            b = s['steps_per_reposition'][mode]
            lines.append(f'\\figRepBox{{{mode}}}{{{x}}}' + ''.join(
                '{' + fmt(b[k]) + '}' for k in ('p5', 'p25', 'p50', 'p75', 'p95', 'mean')))
            tsv.append([tag, run, mode, s['repositions_per_base'], s['repositions']] +
                       [b[k] for k in ('p5', 'p25', 'p50', 'p75', 'p95', 'mean')])
        b = s['ratio_dual_to_min']
        lines.append(f'\\figRepRatio{{{x}}}' + ''.join(
            '{' + fmt(b[k]) + '}' for k in ('p5', 'p25', 'p50', 'p75', 'p95', 'mean')))
        tsv.append([tag, run, 'dual/min', s['repositions_per_base'], s['repositions']] +
                   [b[k] for k in ('p5', 'p25', 'p50', 'p75', 'p95', 'mean')])
    (out/f'{tag}-boxes.tex').write_text('\n'.join(lines) + '\n')
    T = T or tag.upper()
    (out/f'{tag}-ticks.tex').write_text(f'\\def\\FigRepTicks{T}{{{",".join(ticks)}}}\n'
                                        f'\\def\\FigRepTickLabels{T}{{{",".join(labels)}}}\n'
                                        f'\\def\\FigRepXmax{T}{{{len(rows) + 0.6}}}\n')


# Violin layout in batch units: each mode's offset within its batch, and the half width
# at a violin's widest point.
VIOLIN_OFFSET = {'psi': -0.22, 'phi': 0.0, 'dual': 0.22}
VIOLIN_HALF = 0.1
# With --with-min, a fourth violin per batch for the per-repositioning min(psi, phi),
# beside dual, with narrower violins so that four fit.
VIOLIN_OFFSET_MIN = {'psi': -0.3, 'phi': -0.1, 'dual': 0.1, 'min': 0.3}
VIOLIN_HALF_MIN = 0.085


def write_violins(out, tag, runs, violins, with_min=False):
    """Each violin as a closed polygon in axis coordinates, and each mode's means as a line."""
    modes = (*MODES, 'min') if with_min else MODES
    offset = VIOLIN_OFFSET_MIN if with_min else VIOLIN_OFFSET
    half = VIOLIN_HALF_MIN if with_min else VIOLIN_HALF
    lines = [f'% Fig Rep series {tag.upper()} violins: written by bench/campaign-20261012/fig_rep_tikz.py '
             f'from fig_rep_violin.py densities.']
    means = {m: [] for m in modes}
    for x, run in enumerate(runs, start=1):
        for mode in modes:
            v = violins[run][mode]
            cx = x + offset[mode]
            ys = [10 ** t for t in v['log10_steps']]
            right = [f'({cx + half * d:.4f},{y:.5g})' for y, d in zip(ys, v['density'])]
            left = [f'({cx - half * d:.4f},{y:.5g})' for y, d in zip(reversed(ys), reversed(v['density']))]
            lines.append(f'\\figRepViolin{{{mode}}}{{' + ' '.join(right + left) + '}')
            means[mode].append(f'({cx:.4f},{v["mean"]:.5g})')
    for mode in modes:
        lines.append(f'\\figRepMeanLine{{{mode}}}{{' + ' '.join(means[mode]) + '}')
    (out/f'{tag}{"m" if with_min else ""}-violins.tex').write_text('\n'.join(lines) + '\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--a-dir', required=True)
    ap.add_argument('--b-dir', required=True)
    ap.add_argument('--out', default=str(OUT))
    ap.add_argument('--violins', help="fig_rep_violin.py's JSON, to also write the violin files")
    ap.add_argument('--design', choices=sorted(DESIGNS), default='1')
    ap.add_argument('--with-min', action='store_true',
                    help='also write <tag>m-violins.tex with a fourth, min(psi, phi), violin per batch')
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tsv = []
    prefix, series_a, series_b, (ta, tb), (ma, mb) = DESIGNS[args.design]

    def top(s):
        # Designs after the first put repositions per base on the first label line of
        # both series, in the figRate color that fig-rep-v5-d*.tex defines.
        r = rate(s['repositions_per_base'])
        # Leading zeros are dropped so that nine labels fit across series B.
        return r if args.design == '1' else f'\\textcolor{{figRate}}{{{nolead(r)}}}'

    def nolead(t):
        return t[1:] if t.startswith('0.') else t
    a = Path(args.a_dir)
    rows = []
    for run, err in series_a:
        s = json.loads((a/f'runs/{prefix}{run}/summary.json').read_text())
        rows.append((run, f"{top(s)}\\\\{err if args.design == '1' else nolead(err)}", s))
    write(out, ta, rows, tsv, ma)
    b = Path(args.b_dir)
    rows = []
    for cell in series_b:
        path = b/f'runs/chr22-B-{cell}/summary.json'
        if not path.exists():
            print('missing', path)
            continue
        s = json.loads(path.read_text())
        rows.append((cell, f"{top(s)}\\\\{cell[1:]}" if args.design != '1' else
                     f"{cell[1:]}\\\\{rate(s['repositions_per_base'])}", s))
    write(out, tb, rows, tsv, mb)
    if args.violins:
        v = json.loads(Path(args.violins).read_text())['runs']
        write_violins(out, ta, [f'{prefix}{run}' for run, _ in series_a], v)
        write_violins(out, tb, [f'chr22-B-{c}' for c in series_b if f'chr22-B-{c}' in v], v)
        if args.with_min:
            write_violins(out, ta, [f'{prefix}{run}' for run, _ in series_a], v, True)
            write_violins(out, tb, [f'chr22-B-{c}' for c in series_b if f'chr22-B-{c}' in v], v, True)
    with (out/f'batches-{args.design}.tsv' if args.design != '1' else out/'batches.tsv').open('w') as f:
        f.write('series\trun\tmode\trepositions_per_base\trepositions\tp5\tp25\tp50\tp75\tp95\tmean\n')
        for r in tsv:
            f.write('\t'.join(fmt(v) if isinstance(v, float) else str(v) for v in r) + '\n')
    print(out)


if __name__ == '__main__':
    main()
