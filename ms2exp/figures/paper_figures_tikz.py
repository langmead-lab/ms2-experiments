#!/usr/bin/env python3
"""Data step for the TikZ/pgfplots versions of the ms2 paper figures.

Usage:
  paper_figures_tikz.py [--csv FILE] [--captions FILE] [--out DIR]

Reads plotted-points.csv and CAPTIONS.md as written by paper_figures.py (run
that first after new runs; defaults figures/draft/) and writes, under DIR (default
figures/tikz/data):

  <fig>/<panel>.dat         one row per plotted point: key layout x y dx
  <fig>/<panel>-labels.dat  direct labels: x y dx dy anchor text
  <fig>/limits.tex          axis limits shared by the panels of a figure family
  captions.tex              figure captions from CAPTIONS.md, TeX-escaped
  fig-ops/<reference>.dat   Fig Ops bars, one row per layout:
                            key x ns oracle lf total nslab totlab
                            (ns per base for the ordinary query and its
                            oracle replay; LF and total operations per base;
                            preformatted bar labels)
  fig-ops/limits.tex        Fig Ops time axis maximum and tick step per reference
  fig-loracle/, fig-loracle1/
                            Fig LOracle at N = 100 and N = 1, as the locate figures
                            plus psiO and dualO oracle rows at the ordinary point's
                            peak and dodge
  fig-ops2/                 the same for Fig Ops2 (simulated reads from the held-out
                            genomes, of the type paper_figures.py --ops2-reads drew),
                            once plotted-points.csv has its rows;
                            its limits.tex macros begin FigOpsTwo and also give the
                            operations axis maximum and step

The TeX side (figures-style.tex and fig-*.tex) reads only these files, so
regenerating after new runs is: paper_figures.py, then this script, then
latexmk in figures/tikz; python -m ms2exp.figures runs all three.

Columns in <panel>.dat:
  key     series key; figures-style.tex has one pgfplots style per key
          (lcp, psi, dual; lcpO and lcpOS for oracle and oracle-scan;
          lcpF for ftab; rl1, rl19, rl25, rl30, rl40 for ropebwt3 MS-only;
          moni, phoni, augphoni; s11, s8, s6 for ropebwt3 locate)
  layout  rows sharing a layout are joined by a thin line (an ms2 layout's
          oracle ticks or ftab variant, the ropebwt3 -s series)
  x       peak memory in GiB: the plotted value from the CSV for MS figures
          (oracle rows sit at the ordinary point's peak) and the measured
          peak for locate figures, whose dodge is the dx column instead
  y       Mbases/s (MS) or reported locations/s (locate)
  dx      horizontal dodge in pt applied to the marker and its label, so
          points that share a peak stay visible without moving the data
"""
import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
DEFAULT_CSV = ROOT / 'figures' / 'draft' / 'plotted-points.csv'
DEFAULT_CAPTIONS = ROOT / 'figures' / 'draft' / 'CAPTIONS.md'
DEFAULT_OUT = ROOT / 'figures' / 'tikz' / 'data'

MS_FIGS = ['fig-ms', 'fig-ms-oracle', 'fig-ms-ftab']
CHR22_FIGS = ['fig-ms22-10', 'fig-ms22-100']
LOC_FIGS = ['fig-l1', 'fig-l10', 'fig-l100']
LFTAB = 'fig-lftab'
# Fig LOracle: locate at N = 100 and N = 1 with oracle ticks for ms2-psi and ms2-dual.
LORACLE_FIGS = ['fig-loracle', 'fig-loracle1']
OPS = 'fig-ops'
OPS2 = 'fig-ops2'
ALL_FIGS = MS_FIGS + CHR22_FIGS + LOC_FIGS + [LFTAB, OPS]

# Fig Ops: bar order and the TeX macro suffix for each reference's limits.
# Fig Ops time panels show all four layouts; its operation panels leave out ms2-min (figures-style.tex).
# Bar order puts ms2-dual in the rightmost slot.
OPS_KEYS = ['lcp', 'psi', 'phi', 'dual']
OPS_REF_MACRO = {'mtb152': 'Mtb', 'human472': 'Human'}
OPS_STEPS = [10, 20, 25, 50, 100, 200, 250, 500]
# Fig Ops2: limits macro prefix and operations axis steps (read queries reposition far more often).
OPS_MACRO = {OPS: 'FigOps', OPS2: 'FigOpsTwo'}
OPS_COUNT_STEPS = [1, 2, 5, 10, 20, 50, 100, 200]

LAYOUT_KEY = {'ms': 'lcp', 'psi': 'psi', 'dual': 'dual', 'phi': 'phi'}
# The CSV names ms2-min by its mode; the TeX styles key it lcp.
SERIES_KEY = {'min': 'lcp'}
MODE_SUFFIX = {'ordinary': '', 'oracle': 'O', 'oracle-scan': 'OS', 'ftab': 'F'}
ROPE_KEY = {'rope-l1': 'rl1', 'rope-default': 'rl19', 'rope-l25': 'rl25', 'rope-l30': 'rl30',
            'rope-l40': 'rl40'}

# Horizontal dodge (pt) for locate points; ms2-psi, ms2-phi and ms2-dual share
# one index file and so one peak.
LOC_DODGE = {'psi': -3.0, 'phi': 0.0, 'dual': 3.0}
LFTAB_DODGE = {'psi': -14.0, 'psiF': -9.0, 'phi': -2.5, 'phiF': 2.5, 'dual': 9.0, 'dualF': 14.0}
# Fig LOracle spreads psi and dual wider so that their oracle ticks stay apart; each oracle row
# shares its ordinary point's peak and dodge, so the tick sits on the vertical line above the point.
LORACLE_DODGE = {'psi': -3.5, 'psiO': -3.5, 'dual': 3.5, 'dualO': 3.5}
# ropebwt3 MS-only modes share one peak (within about 1 GiB) and differ in
# speed by as little as 4 percent, so they fan out to the right of -l1.
ROPE_FAN = {'rl1': 0.0, 'rl19': 2.7, 'rl25': 5.4, 'rl30': 8.1, 'rl40': 10.8}
FAN_NOTE = ('ropebwt3 modes other than -l1 are shifted right in steps of 2.7 pt (a few GiB) so that '
            'points sharing one peak memory stay visible; plotted-points.csv keeps the measured peaks.')
LOC_NOTE = ('ms2-psi, ms2-phi and ms2-dual are drawn a few points apart around their shared peak; '
            'plotted-points.csv keeps the measured peaks.')

# Direct labels: (key, text, anchor, dx, dy, panels). panels is 'first'
# (leftmost panel only) or 'all'. Offsets are in pt, added to the point's dodge.
MS_LABELS = [('lcp', 'min', 'east', -2.2, 0, 'first'), ('psi', 'psi', 'east', -2.2, 0, 'first'),
             ('dual', 'dual', 'east', -2.2, 0, 'first'), ('rl1', '-l1', 'north west', -2.0, -2.0, 'first')]
# The chr22 ms2-phi label is drawn with a leader line by gap_annotations.py, clear of the MONI cluster.
CHR22_LABELS = MS_LABELS + [('moni', 'MONI', 'south', 0, 2.2, 'first'),
                            ('phoni', 'PHONI', 'north', 0, -2.2, 'first'),
                            ('augphoni', 'AugPHONI', 'west', 2.5, 0, 'first')]
LOC_LABELS = [('psi', 'psi, dual', 'east', -2.5, 0, 'first'), ('s11', '-s11', 'west', 3.0, 0, 'first'),
              ('s8', '-s8', 'north', 0, -2.2, 'first'), ('s6', '-s6', 'north', 0, -2.2, 'first')]
# At N = 1 the ropebwt3 points lie nearly flat, so -s8 and -s6 go above their points.
LOC_LABELS_FLAT = LOC_LABELS[:2] + [('s8', '-s8', 'south', 0, 2.2, 'first'), ('s6', '-s6', 'south', 0, 2.2, 'first')]
# Fig L1 also carries a speed-gap arrow rising from -s6, so its -s6 label sits below and left.
LOC_LABELS_L1 = LOC_LABELS_FLAT[:3] + [('s6', '-s6', 'north east', -1.5, -2.2, 'first')]
FIG_LABELS = {'fig-ms': MS_LABELS, 'fig-ms-oracle': MS_LABELS, 'fig-ms-ftab': MS_LABELS,
              'fig-ms22-10': CHR22_LABELS, 'fig-ms22-100': CHR22_LABELS,
              'fig-l1': LOC_LABELS_L1, 'fig-l10': LOC_LABELS, 'fig-l100': LOC_LABELS, LFTAB: [],
              'fig-loracle': LOC_LABELS, 'fig-loracle1': LOC_LABELS_FLAT}
# fig-lftab names its points in the legend and labels only each ftab point, with its speed ratio over the
# same layout without ftab.
RATIO_LABELS = {'psiF': ('south east', 0, 2.0), 'phiF': ('north east', -3.0, -1.5), 'dualF': ('south west', 0, 2.0)}

PANEL_ORDER = ['human_short', 'human_pacbio', 'human_ont_modern', 'guccione_skin',
               'simshort', 'simpacbio', 'simont', 'cap1', 'cap10', 'cap100']


def series_key(cell):
    """Maps a summary cell id to (key, layout, cap)."""
    if cell in ROPE_KEY:
        return ROPE_KEY[cell], ROPE_KEY[cell], None
    if cell in ('moni', 'phoni', 'augphoni'):
        return cell, cell, None
    m = re.fullmatch(r'cap(\d+)-rope-s(\d+)', cell)
    if m:
        return f's{m.group(2)}', 'rs', f'cap{m.group(1)}'
    m = re.fullmatch(r'(?:cap(\d+)-)?(ms|psi|dual|phi)-(ordinary|oracle-scan|oracle|ftab)', cell)
    if m:
        layout = LAYOUT_KEY[m.group(2)]
        return layout + MODE_SUFFIX[m.group(3)], layout, (f'cap{m.group(1)}' if m.group(1) else None)
    raise ValueError(f'unknown cell id {cell}')


def read_points(csv_path):
    """Returns {fig: {panel: [point dict]}} with keys assigned."""
    figs = defaultdict(lambda: defaultdict(list))
    with open(csv_path) as fh:
        for r in csv.DictReader(fh):
            fig = r['figure']
            if fig in (OPS, OPS2):
                continue
            # Oracle-scan replay is not part of the paper figures.
            if r['cell'].endswith('-oracle-scan'):
                continue
            key, layout, cap = series_key(r['cell'])
            locate = fig in LOC_FIGS or fig == LFTAB or fig in LORACLE_FIGS
            panel = cap if fig == LFTAB else r['read_set']
            x = float(r['peak_gib_measured'] if locate else r['peak_gib_plotted'])
            if fig == LFTAB:
                dx = LFTAB_DODGE.get(key, 0.0)
            elif fig in LORACLE_FIGS:
                dx = LORACLE_DODGE.get(key, 0.0)
            elif locate:
                dx = LOC_DODGE.get(key, 0.0)
            else:
                dx = ROPE_FAN.get(key, 0.0)
            figs[fig][panel].append({'key': key, 'layout': layout, 'x': x, 'y': float(r['speed']), 'dx': dx})
    for fig in LORACLE_FIGS:
        for pts in figs.get(fig, {}).values():
            ordinary = {p['key']: p['x'] for p in pts}
            for p in pts:
                if p['key'].endswith('O') and p['layout'] in ordinary:
                    p['x'] = ordinary[p['layout']]
    return figs


def ordered_panels(panels):
    return sorted(panels, key=PANEL_ORDER.index)


def write_panel(path, pts):
    # Rows of one layout are contiguous and sorted so joining lines run monotonically.
    pts = sorted(pts, key=lambda p: (p['layout'], p['x'], p['y']))
    with open(path, 'w') as fh:
        fh.write('key layout x y dx\n')
        for p in pts:
            fh.write(f"{p['key']} {p['layout']} {p['x']:.4f} {p['y']:.6g} {p['dx']:g}\n")


def write_labels(path, fig, pts, first):
    by_key = {p['key']: p for p in pts}
    rows = []
    for key, text, anchor, dx, dy, panels in FIG_LABELS[fig]:
        if key in by_key and (panels == 'all' or first):
            p = by_key[key]
            rows.append((p['x'], p['y'], p['dx'] + dx, dy, anchor, text))
    if fig == LFTAB:
        for key, (anchor, dx, dy) in RATIO_LABELS.items():
            base = by_key.get(key[:-1])
            if key in by_key and base:
                p = by_key[key]
                rows.append((p['x'], p['y'], p['dx'] + dx, dy, anchor, f"$\\times${p['y'] / base['y']:.2f}"))
    with open(path, 'w') as fh:
        fh.write('x y dx dy anchor text\n')
        for x, y, dx, dy, anchor, text in rows:
            fh.write(f'{x:.4f} {y:.6g} {dx:g} {dy:g} {{{anchor}}} {{{text}}}\n')


def limits(points, ylo, yhi, xhi=1.08):
    xs = [p['x'] for p in points]
    ys = [p['y'] for p in points]
    return max(xs) * xhi, min(ys) / ylo, max(ys) * yhi


def write_limits(path, lim, panels):
    xmax, ymin, ymax = lim
    with open(path, 'w') as fh:
        fh.write('% Written by paper_figures_tikz.py; axis limits shared by this figure family.\n')
        fh.write(f'\\def\\FigXmax{{{xmax:.4g}}}\n\\def\\FigYmin{{{ymin:.4g}}}\n\\def\\FigYmax{{{ymax:.4g}}}\n')
        fh.write(f"\\def\\FigPanels{{{','.join(panels)}}}\n")


def family_points(figs, names, panel_filter=None):
    return [p for f in names if f in figs for panel, pts in figs[f].items()
            if panel_filter is None or panel in panel_filter for p in pts]


def read_ops(csv_path, fig=OPS):
    """Returns {reference: {key: row}} for Fig Ops (or Fig Ops2) from plotted-points.csv."""
    refs = defaultdict(lambda: defaultdict(dict))
    with open(csv_path) as fh:
        for r in csv.DictReader(fh):
            if r['figure'] != fig:
                continue
            m = re.fullmatch(r'ms2-(\w+)( oracle-scan| oracle| ops)?', r['series'])
            key, kind = SERIES_KEY.get(m.group(1), m.group(1)), (m.group(2) or '').strip()
            row = refs[r['read_set']][key]
            if kind == '':
                row['ns'] = float(r['speed'])
            elif kind == 'oracle':
                row['oracle'] = float(r['speed'])
            elif kind == 'ops':
                row['total'] = float(r['speed'])
                row['lf'] = float(re.match(r'LF ([0-9.]+)', r['x_note']).group(1))
    return refs


def nice_step(top):
    return next((s for s in OPS_STEPS if top / s <= 6), OPS_STEPS[-1])


def write_ops(refs, d, fig=OPS):
    d.mkdir(parents=True, exist_ok=True)
    title = 'Fig Ops2' if fig == OPS2 else 'Fig Ops'
    lim = [f'% Written by paper_figures_tikz.py; {title} time axis per reference.']
    prefix = OPS_MACRO[fig]
    for ref, rows in refs.items():
        # A capped cell (Fig Ops2 ms2-phi on human reads) has no time or count row, so it has no bar.
        rows = {k: r for k, r in rows.items() if 'ns' in r and 'total' in r}
        with open(d / f'{ref}.dat', 'w') as fh:
            fh.write('key x ns oracle lf total nslab totlab\n')
            for i, key in enumerate(OPS_KEYS):
                if key not in rows:
                    continue
                r = rows[key]
                fh.write(f"{key} {i} {r['ns']:.3f} {r.get('oracle', float('nan')):.3f} {r['lf']:.5f} "
                         f"{r['total']:.5f} {r['ns']:.0f} {r['total']:.2f}\n")
        top = (max((rows[k]['ns'] for k in OPS_KEYS if k in rows), default=0.0) * 1.16) or 1.0
        step = nice_step(top)
        name = OPS_REF_MACRO[ref]
        lim.append(f'\\def\\{prefix}TimeMax{name}{{{top:.4g}}}\\def\\{prefix}TimeStep{name}{{{step}}}')
    if fig == OPS2:
        top = max((r['total'] for rows in refs.values() for r in rows.values() if 'total' in r), default=4.0)
        step = next((s for s in OPS_COUNT_STEPS if top / s <= 4), OPS_COUNT_STEPS[-1])
        lim.append(f'\\def\\{prefix}CountMax{{{step * -(-top // step):g}}}\\def\\{prefix}CountStep{{{step}}}')
    (d / 'limits.tex').write_text('\n'.join(lim) + '\n')


def tex_escape(s):
    s = s.replace('\\', r'\textbackslash{}')
    for a, b in [('&', r'\&'), ('%', r'\%'), ('#', r'\#'), ('_', r'\_'), ('$', r'\$'), ('~', r'\textasciitilde{}')]:
        s = s.replace(a, b)
    s = s.replace('×', r'$\times$')
    s = re.sub(r'`([^`]*)`', r'\\texttt{\1}', s)
    # Command-line options such as -l1 and -s11 go in typewriter type, where l and 1 differ.
    s = re.sub(r'(?<![\w-])(-[a-zA-Z]\w*)', r'\\texttt{\1}', s)
    s = s.replace(' (Load times table below)', '')
    return s


def adapt_caption(fig, text):
    """Replaces the matplotlib draft's offset wording with this version's dodges."""
    text = re.sub(r'ms2-psi, ms2-phi and ms2-dual are offset horizontally by [0-9.]+ GiB steps around their '
                  r'shared peak for visibility\.', LOC_NOTE, text)
    if fig in MS_FIGS + CHR22_FIGS:
        text = text.replace(' Reference:', ' ' + FAN_NOTE + ' Reference:', 1)
    return text


def write_captions(md_path, out_path):
    text = Path(md_path).read_text()
    caps = {}
    for m in re.finditer(r'^## .*?\(`([^`]+)`\)\n\n!\[[^\]]*\]\([^)]*\)\n\n(.+?)\n', text, re.M):
        caps[m.group(1)] = m.group(2)
    with open(out_path, 'w') as fh:
        fh.write('% Written by paper_figures_tikz.py from CAPTIONS.md. Use \\FigCaption{fig-ms}.\n')
        fh.write('\\newcommand{\\FigCaption}[1]{\\csname figcaption:#1\\endcsname}\n')
        # Fig Ops2 and Fig LOracle get captions once CAPTIONS.md has them.
        for fig in ALL_FIGS + [f for f in [OPS2] + LORACLE_FIGS if f in caps]:
            body = adapt_caption(fig, caps.get(fig, 'Caption not found in CAPTIONS.md.'))
            body = tex_escape(body)
            fh.write(f'\\expandafter\\def\\csname figcaption:{fig}\\endcsname{{{body}}}\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--csv', default=DEFAULT_CSV)
    ap.add_argument('--captions', default=DEFAULT_CAPTIONS)
    ap.add_argument('--out', default=DEFAULT_OUT)
    args = ap.parse_args()
    out = Path(args.out)
    figs = read_points(args.csv)
    lims = {}
    lims.update({f: limits(family_points(figs, MS_FIGS), 1.6, 1.8) for f in MS_FIGS})
    # Fig MS+ftab shows only the short reads, so it gets its own vertical range rather than the
    # family's, whose floor is set by the slowest long-read points.
    lims['fig-ms-ftab'] = limits(family_points(figs, ['fig-ms-ftab']), 1.6, 1.8)
    for f in CHR22_FIGS:
        # A deeper floor keeps the slowest ms2-phi point clear of the panel titles.
        lims[f] = limits(family_points(figs, CHR22_FIGS), 2.8, 1.8)
    lims.update({f: limits(family_points(figs, LOC_FIGS + [LFTAB]), 3.2, 2.0) for f in LOC_FIGS})
    # Figs L1 and L10 report far fewer locations per second than the rest of the family, so each
    # gets its own vertical range rather than leaving its upper decades empty.
    for f in ['fig-l1', 'fig-l10']:
        lims[f] = limits(family_points(figs, [f]), 3.2, 2.0)
    # Extra room on the right keeps the six spread ms2 points and their ratio labels inside the panel.
    lims[LFTAB] = limits(family_points(figs, [LFTAB]), 3.2, 3.2, xhi=1.32)
    for f in LORACLE_FIGS:
        lims[f] = limits(family_points(figs, LOC_FIGS + [LFTAB] + LORACLE_FIGS), 3.2, 2.0)
    for fig, panels in figs.items():
        d = out / fig
        d.mkdir(parents=True, exist_ok=True)
        order = ordered_panels(panels)
        for i, panel in enumerate(order):
            write_panel(d / f'{panel}.dat', panels[panel])
            write_labels(d / f'{panel}-labels.dat', fig, panels[panel], first=(i == 0))
        write_limits(d / 'limits.tex', lims[fig], order)
    ops = read_ops(args.csv)
    if ops:
        write_ops(ops, out / OPS)
    ops2 = read_ops(args.csv, OPS2)
    if ops2:
        write_ops(ops2, out / OPS2, OPS2)
    write_captions(args.captions, out / 'captions.tex')
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
