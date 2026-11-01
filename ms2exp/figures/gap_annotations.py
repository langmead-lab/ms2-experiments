#!/usr/bin/env python3
"""Speed-gap annotations for the speed-against-memory figures in figures-tikz.

For each panel, a dashed line runs at the speed of the fastest ordinary ms2 point, from that
point to the x position of a ropebwt3 reference point, and an arrow rises from the reference
point to the line, labeled with the speed ratio. The reference is ropebwt3 -l1 for the MS
figures and ropebwt3 -s6 (its fastest locate mode) for Fig L100. Writes one TikZ snippet per
panel to DATA/gap/<figure>-<panel>.tex, drawn inside that panel's axis. Run it
after paper_figures_tikz.py, since it reads that script's panel data.
"""
import argparse
from pathlib import Path

# figure: (panels, reference key, ms2 keys, label side of the arrow, label position along it)
FIGS = {
    'fig-ms22-100': (['simshort', 'simpacbio', 'simont'], 'rl1', ['lcp', 'psi', 'dual'], 'east', 0.8),
    'fig-ms22-10': (['simshort', 'simpacbio', 'simont'], 'rl1', ['lcp', 'psi', 'dual'], 'east', 0.8),
    'fig-ms-oracle': (['human_short', 'human_pacbio', 'human_ont_modern', 'guccione_skin'], 'rl1',
                      ['lcp', 'psi', 'dual'], 'east', 0.8),
    'fig-l1': (['human_short', 'human_pacbio', 'human_ont_modern', 'guccione_skin'], 's6',
               ['psi', 'phi', 'dual'], 'west', 0.5),
    'fig-l10': (['human_short', 'human_pacbio', 'human_ont_modern', 'guccione_skin'], 's6',
                ['psi', 'phi', 'dual'], 'west', 0.5),
    'fig-l100': (['human_short', 'human_pacbio', 'human_ont_modern', 'guccione_skin'], 's6',
                 ['psi', 'phi', 'dual'], 'west', 0.5),
}
# Figures whose first panel draws no point labels and leaves the naming to the legend.
NO_POINT_LABELS = {'fig-ms22-100', 'fig-ms22-10', 'fig-ms-oracle', 'fig-l1', 'fig-l10', 'fig-l100'}
# Panels whose ratio label would cover an ms2 point beside the arrow: (label anchor, position along
# the arrow). There the label sits above the arrow tip instead.
LABEL_AT = {('fig-l1', 'human_short'): ('south east', 1.0), ('fig-l10', 'human_short'): ('south east', 1.0)}


def points(path):
    rows = path.read_text().split('\n')[1:]
    return {f[0]: dict(x=float(f[2]), y=float(f[3]), dx=float(f[4])) for f in (r.split() for r in rows) if f}


def ratio_text(r):
    return f'{r:.0f}' if r >= 10 else f'{r:.1f}'


def write(data):
    """Write DATA/gap/<figure>-<panel>.tex from the panel files paper_figures_tikz.py wrote under DATA."""
    DATA = Path(data)
    out = DATA / 'gap'
    out.mkdir(exist_ok=True)
    for fig, (panels, ref, keys, side, pos) in FIGS.items():
        for panel in panels:
            p = points(DATA / fig / f'{panel}.dat')
            r, m = p[ref], max((p[k] for k in keys if k in p), key=lambda q: q['y'])
            # The label sits on the side of the arrow away from the ms2 points.
            label_anchor = 'west' if side == 'east' else 'east'
            label_anchor, label_pos = LABEL_AT.get((fig, panel), (label_anchor, pos))
            # The first panel carries the point labels; in the MS figures the ms2-min label sits
            # left of its point, on the dashed line, so the line stops short of it there.
            trim = ', shorten <=15pt' if panel == panels[0] and side == 'east' and fig not in NO_POINT_LABELS else ''
            (out / f'{fig}-{panel}.tex').write_text(
                f'% Written by gap_annotations.py.\n'
                f'\\draw[fig gap line{trim}] (axis cs:{m["x"]},{m["y"]}) -- (axis cs:{r["x"]},{m["y"]});\n'
                f'\\draw[fig gap arrow] (axis cs:{r["x"]},{r["y"]}) -- (axis cs:{r["x"]},{m["y"]})\n'
                f'  node[fig gap label, anchor={label_anchor}, pos={label_pos}] {{\\texttimes{{}}{ratio_text(m["y"] / r["y"])}}};\n')
    print(f'wrote {out}')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('data', help="paper_figures_tikz.py's output directory")
    write(ap.parse_args().data)


if __name__ == '__main__':
    main()
