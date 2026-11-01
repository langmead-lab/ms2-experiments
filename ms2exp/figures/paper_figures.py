#!/usr/bin/env python3
"""Draft paper figures for the ms2 (Orbit) matching-statistics paper.

Usage:
  paper_figures.py [--ms-dir DIR] [--chr22-dir DIR] [--locate-dir DIR] [--ops-dir DIR]
                   [--ops2-reads pacbio|ont] [--out DIR]

Reads the exploratory Human472 and chr22 summaries (explore_ms.py,
explore_ms_chr22.py, explore_offsets.py and explore_ops.py summary.json files) and writes, as
PNG, SVG and PDF on a white background:

  fig-ms          MS-only speed against peak memory, ordinary ms2 layouts and
                  every ropebwt3 MS-only mode, one panel per read set
  fig-ms-oracle   fig-ms plus oracle replay for each ms2 layout
  fig-ms-ftab     short reads only, ms2 layouts with and without ftab
  fig-ms22-10, fig-ms22-100
                  MS-only on chr22-h10 and chr22-h100 with simulated short,
                  PacBio HiFi and Nanopore reads, adding MONI, PHONI, AugPHONI
  fig-l1, fig-l10, fig-l100
                  SMEM offset reporting (locate) at caps 1, 10 and 100
  fig-lftab       short reads only, the three caps with and without ftab
  fig-loracle, fig-loracle1
                  locate at N = 100 and N = 1 with oracle replay ticks for
                  ms2-psi and ms2-dual
  fig-ops         time and operations per base for ms2-lcp, ms2-psi, ms2-dual
                  and ms2-phi on held-out genomes (mtb152 and HPRC2 full), with
                  oracle replay
  fig-ops2        as fig-ops on simulated reads drawn from the same held-out
                  genomes (chr22 recipe; --ops2-reads picks PacBio HiFi, the
                  default, or Nanopore); drawn once its summaries are in --ops-dir

It also writes CAPTIONS.md (draft captions with workload sizes taken from the
summaries, a load-time table, and the gaps found in the data) and
plotted-points.csv (every plotted value with its load time and node).

The chr22 summaries are copied from Rockfish into --chr22-dir as
chr22-h{10,100}-ms-sim{short,pacbio,ont}/summary.json (source runs:
results/campaign-20261012/explore/chr22-h*-ms-sim*-e637ac0-v1). A chr22 figure
draws only the panels whose summary is present; copy a missing summary in and
rerun to complete it.

Conventions shared by every figure: one color per ms2 layout (lcp, psi, dual),
one marker per mode (ordinary, ftab). The ms2 MS-only set is ms2-lcp, ms2-psi
and ms2-dual. ropebwt3 -l1 and the ropebwt3 SMEM reporting set (-l1 at SSA
shifts 11, 8 and 6) are red. ropebwt3 at longer minimum lengths reports only
SMEMs of at least that length; those points are a gray ramp and their labels
carry an asterisk. MONI, PHONI and AugPHONI share one purple and are hollow,
because their speed comes from the one-read wall fallback rather than the
load/query split. Speed is on a log axis and measures query time only; load
time is excluded. Oracle replay is a horizontal tick on a vertical line rising
from the layout's ordinary point, at the ordinary point's peak memory, because
replay loads only LF and its own peak would understate what a real query needs.
Runs may also record oracle-scan replay; no figure, caption or CSV row uses it.
"""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
CAMPAIGN = ROOT / 'reports' / 'campaign-20261012'
DEFAULT_MS = CAMPAIGN / 'explore-ms' / 'summaries-v2'
DEFAULT_CHR22 = CAMPAIGN / 'explore-ms' / 'summaries-v2'
DEFAULT_LOCATE = CAMPAIGN / 'explore-oracle-locate' / 'summaries'
DEFAULT_OPS = CAMPAIGN / 'explore-ops' / 'summaries'
DEFAULT_OUT = CAMPAIGN / 'figures-draft'

COHORTS = ['human_short', 'human_pacbio', 'human_ont_modern', 'guccione_skin']
CHR22_REFS = [10, 100]
CHR22_SETS = ['simshort', 'simpacbio', 'simont']
CHR22_RUN = 'chr22-h{h}-ms-{s}'
CHR22_REMOTE = 'results/campaign-20261012/explore/chr22-h{h}-ms-{s}-e637ac0-v1/summary.json'
TITLES = {
    'human_short': 'Short reads (Illumina, ACGT-only)',
    'human_pacbio': 'PacBio HiFi',
    'human_ont_modern': 'Nanopore (HG002 Q28)',
    'guccione_skin': 'Skin metagenome (Illumina, ERR14066769)',
    'simshort': 'Simulated short reads',
    'simpacbio': 'Simulated PacBio HiFi',
    'simont': 'Simulated Nanopore',
    'mtb152': 'mtb152, 8 held-out genomes and their reverse complements',
    'human472': 'HPRC2 full, APR008 haplotype 1',
}
# Reference panels, from data/ladder/chr22-h*/manifest.json as recorded in the
# held-out haplotype manifest.
CHR22_REF_TEXT = {10: 'chr22-h10 (chr22 of 10 HPRC haplotypes from 5 donors)',
                  100: 'chr22-h100 (chr22 of 100 HPRC haplotypes from 50 donors)'}
# Simulation details, from inputs/chr22_sim/{short,pacbio,ont}-60mb-v1/manifest.json.
CHR22_SIM_TEXT = {
    'simshort': ('Simulated short reads: ART 2.5.8, HiSeq 2500 profile (-ss HS25), single-end 125 bp'),
    'simpacbio': ('Simulated PacBio HiFi: pbsim3 3.0.5 multi-pass subreads (QSHMM-RSII model, 10 passes) made '
                  'into HiFi reads by ccs 6.4.0 (defaults: at least 3 passes, predicted accuracy at least 0.99), '
                  'with read lengths drawn from the real PacBio HiFi set (median 20,177 bases)'),
    'simont': ('Simulated Nanopore: pbsim3 3.0.5 with its standard ONT model (ERRHMM-ONT-HQ) at mean accuracy '
               '0.96, the model\'s ceiling (measured 95.98 percent), below the accuracy of real HG002 Q28 data, '
               'with read lengths drawn from the real Nanopore (HG002 Q28) set (N50 94,132 bases)'),
}
# Read counts from the same manifests, used only when a summary is not yet present.
CHR22_SIM_SIZE = {'simshort': (480000, 60000000), 'simpacbio': (2867, 60005589), 'simont': (2682, 60033032)}
MS_RUN = '{cohort}-ms-acgt60mb'
# Locate runs: the October 25 oracle-locate runs (explore_offsets.py prepare --oracle), which add oracle
# replay to the October 23 s6 cells and run ropebwt3 with -K 1m.
LOCATE_RUNS = {
    'human_short': 'human_short-acgt49152-e637ac0-oracle-v2',
    'human_pacbio': 'human_pacbio-acgt20mb-e637ac0-oracle-v2',
    'human_ont_modern': 'human_ont_modern-acgt10mb-e637ac0-oracle-v1',
    'guccione_skin': 'guccione_skin-acgt400kb-e637ac0-oracle-v2',
}
# Fig LOracle: locate with oracle replay, at these caps (stem fig-loracle for N = 100, fig-loracle1 for N = 1).
LORACLE_CAPS = {100: 'fig-loracle', 1: 'fig-loracle1'}
CAPS = [1, 10, 100]
# Fig Ops runs (results/campaign-20261012/explore/<run>/summary.json on Rockfish).
OPS_REFS = ['mtb152', 'human472']
OPS_RUNS = {'mtb152': 'fig_ops-mtb152-e637ac0-v1', 'human472': 'fig_ops-human472-apr008-e637ac0-v2'}
OPS_TITLE = {'mtb152': 'mtb152', 'human472': 'HPRC2 full'}
# Fig Ops2 runs: the same cells on simulated reads from the same held-out genomes, one set of runs per read
# type. The Nanopore runs were stopped after mtb152 (EXECUTION.md, 2026-10-26); PacBio HiFi is the default.
OPS2_RUNS = {
    'pacbio': {'mtb152': 'fig_ops2-mtb152-pacbio-e637ac0-v1', 'human472': 'fig_ops2-human472-apr008-pacbio-e637ac0-v1'},
    'ont': {'mtb152': 'fig_ops2-mtb152-ont-e637ac0-v1', 'human472': 'fig_ops2-human472-apr008-ont-e637ac0-v1'},
}
OPS2_DEFAULT_READS = 'pacbio'
# How each read type was simulated, from inputs/fig_ops2/sim-*/manifest.json (the chr22 recipes).
OPS2_SIM_TEXT = {
    'pacbio': ('simulated PacBio HiFi reads in place of whole held-out genomes. The reads are simulated as for '
               'chr22, from the held-out genomes: pbsim3 3.0.5 multi-pass subreads (QSHMM-RSII model, 10 passes) '
               'made into HiFi reads by ccs 6.4.0 (defaults: at least 3 passes, predicted accuracy at least 0.99), '
               'read lengths drawn from the real PacBio HiFi set, both strands, whole ACGT-only reads to about '
               '60 Mbases.'),
    'ont': ('simulated Nanopore reads in place of whole held-out genomes. The reads are simulated as for chr22, '
            'from the held-out genomes: pbsim3 3.0.5 with its standard ONT model (ERRHMM-ONT-HQ) at mean accuracy '
            '0.96, read lengths drawn from the real Nanopore (HG002 Q28) set, both strands, whole ACGT-only reads '
            'to about 60 Mbases.'),
}
HUMAN472 = 'HPRC2 full (Human472, 472 haplotypes)'

# Ink and surface. Text never wears a series color.
INK = '#0b0b0b'
INK_2 = '#52514e'
INK_MUTED = '#8a8984'
GRID = '#e4e3df'
SURFACE = '#ffffff'

# Layout colors, validated as a set with ropebwt3 red and the Others purple
# (validate_palette.js, light, all pairs): no normal-vision failure; the CVD
# warn band is covered by markers, direct labels and the legend.
LAYOUTS = ['ms', 'psi', 'dual']
LAYOUT_NAME = {'ms': 'ms2-lcp', 'psi': 'ms2-psi', 'dual': 'ms2-dual', 'phi': 'ms2-phi'}
# The chr22 figures also draw ms2-phi, which no other speed-against-memory figure shows.
CHR22_PHI_REFS = (10, 100)
LAYOUT_COLOR = {'ms': '#2a78d6', 'psi': '#1baf7a', 'dual': '#c98500', 'phi': '#a3478f'}
ROPE_RED = '#e34948'
# ms2-phi appears only in Fig Ops, so the layout set there adds a plum. The
# four layout colors pass validate_palette.js as a set (light, all pairs). The
# plum sits too close to the Others purple and ropebwt3 red to share a figure
# with them; no figure draws those together.
LAYOUT_COLOR_OPS = dict(LAYOUT_COLOR, phi='#a3478f')
OTHERS_COLOR = '#7d3c98'

MODES = ['ordinary', 'ftab']
MODE_MARKER = {'ordinary': 'o', 'ftab': 'D'}
MODE_SUFFIX = {'ordinary': '', 'oracle': ' oracle', 'ftab': '+ftab'}
# Replay tick (marker size and edge width in points). Only the full oracle is drawn.
REPLAY_MODES = ['oracle']
REPLAY_TICK = {'oracle': (9.0, 1.7)}

# ropebwt3 MS-only modes: -l1 finds every MS value; longer minimum lengths do not.
ROPE_MS = [('rope-l1', 1), ('rope-default', 19), ('rope-l25', 25), ('rope-l30', 30), ('rope-l40', 40)]
ROPE_LONG_STYLE = {19: ('#4d4d4d', 'X'), 25: ('#6e6e6e', 'v'), 30: ('#919191', 'p'), 40: ('#b3b3b3', 'h')}
# ropebwt3 SMEM reporting set: -l1 -pN at SSA shifts 11, 8 and 6.
ROPE_SHIFTS = [11, 8, 6]
ROPE_SHIFT_MARKER = {11: 'P', 8: '*', 6: 'H'}
# Others set: (cell id, name, marker, label offset dx, dy, alignment).
OTHERS = [('moni', 'MONI', 's', 0, 7, 'center'), ('phoni', 'PHONI', 'D', 0, -7, 'center'),
          ('augphoni', 'AugPHONI', '^', 6, 0, 'left')]

MARKER_SIZE = 34
FIG_W = 7.0
ASTERISK_NOTE = ('An asterisk marks ropebwt3 modes that report only SMEMs of at least the given length '
                 '(-l19 is the default), so they cannot find all MS values or all SMEMs.')
ORACLE_NOTE = ('Oracle replay runs a recorded query and loads only LF, so it is drawn at the peak memory of '
               'the same layout\'s ordinary run.')
OTHERS_NOTE = ('MONI, PHONI and AugPHONI (purple) are drawn hollow because their speed uses the one-read wall '
               'fallback: whole-process wall time minus the wall time of a one-read run, rather than the '
               'load/query split that ms2 and ropebwt3 report, so their speed includes any load work that '
               'the one-read subtraction does not cancel.')

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 7.5,
    'axes.titlesize': 8,
    'axes.labelsize': 7.5,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'axes.edgecolor': INK_MUTED,
    'axes.labelcolor': INK_2,
    'axes.titlecolor': INK,
    'xtick.color': INK_2,
    'ytick.color': INK_2,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
})


# Loading

def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def load_ms(ms_dir):
    out = {}
    for cohort in COHORTS:
        path = Path(ms_dir) / MS_RUN.format(cohort=cohort) / 'summary.json'
        if path.exists():
            out[cohort] = load_json(path)
    return out


def load_chr22(chr22_dir):
    """Returns {(h, set): summary} for the chr22 simulated-read runs present."""
    out = {}
    for h in CHR22_REFS:
        for s in CHR22_SETS:
            path = Path(chr22_dir) / CHR22_RUN.format(h=h, s=s) / 'summary.json'
            if path.exists():
                out[(h, s)] = load_json(path)
    return out


def load_locate(locate_dir):
    out = {}
    for cohort, run in LOCATE_RUNS.items():
        path = Path(locate_dir) / run / 'summary.json'
        if path.exists():
            out[cohort] = load_json(path)
    return out


def load_ops(ops_dir, runs=OPS_RUNS):
    out = {}
    for ref, run in runs.items():
        path = Path(ops_dir) / run / 'summary.json'
        if path.exists():
            out[ref] = load_json(path)
    return out


def ops_cells(summary):
    """Fig Ops cells keyed layout-mode; the runs use one interleave (32), so the -k32 suffix is dropped."""
    return {f"{c['layout']}-{c['mode']}": c for c in summary['cells'] if c.get('interleave') == 32}


def cells(summary):
    return {c['id']: c for c in summary['cells']}


def ok(cell):
    return cell is not None and cell.get('status') == 'checked'


def gib(cell):
    return cell['peak_rss_kib'] / 2 ** 20


def workload(summary):
    w = summary['config']['workload']
    return int(w['reads']), int(w['bases'])


def revisions(summaries):
    return sorted({s['config']['orbit']['revision'][:7] for s in summaries})


# The Human472 MS grids (fig-ms, fig-ms-oracle) draw ms2-phi beside ms2-psi.
MS_GRID_LAYOUTS = ['ms', 'psi', 'phi', 'dual']
# Cell ids that set the axis range of each figure family, so fig-ms and
# fig-ms-oracle share axes and out-of-scope cells never stretch them.
HUMAN_MS_IDS = ([f'{l}-{m}' for l in MS_GRID_LAYOUTS for m in ('ordinary', 'oracle', 'ftab')] +
                [cid for cid, _ in ROPE_MS])
CHR22_MS_IDS = [f'{l}-ordinary' for l in LAYOUTS] + [cid for cid, _ in ROPE_MS] + [o[0] for o in OTHERS]


# Drawing helpers

class Recorder:
    """Collects every plotted value so it can be written to plotted-points.csv."""

    def __init__(self):
        self.rows = []

    def add(self, figure, cohort, label, x, y, metric, cell, x_note=''):
        self.rows.append({
            'figure': figure, 'read_set': cohort, 'series': label, 'peak_gib_plotted': round(x, 3),
            'peak_gib_measured': round(gib(cell), 3), 'speed': y, 'metric': metric,
            'median_query_s': cell.get('median_query_seconds'), 'median_load_s': cell.get('median_load_seconds'),
            'x_note': x_note, 'cell': cell['id'],
        })


def style_axes(ax, ylabel, xlabel):
    ax.set_yscale('log')
    ax.grid(True, which='major', color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.yaxis.set_minor_formatter(NullFormatter())


def decade_formatter():
    def fmt(v, _):
        if v >= 1e6:
            return f'{v / 1e6:g}M'
        if v >= 1e3:
            return f'{v / 1e3:g}k'
        return f'{v:g}'
    return FuncFormatter(fmt)


def point(ax, x, y, color, marker, hollow=False, size=MARKER_SIZE, z=3):
    if hollow:
        ax.scatter([x], [y], s=size, marker=marker, facecolors=SURFACE, edgecolors=color, linewidths=1.3, zorder=z)
    else:
        ax.scatter([x], [y], s=size, marker=marker, color=color, edgecolors=SURFACE, linewidths=0.7, zorder=z)


def replay_tick(ax, x, y, color, mode):
    size, width = REPLAY_TICK[mode]
    ax.plot([x], [y], linestyle='none', marker='_', markersize=size, markeredgewidth=width, color=color, zorder=4)


def direct_label(ax, x, y, text, dx=5, dy=0, ha='left'):
    ax.annotate(text, (x, y), xytext=(dx, dy), textcoords='offset points', ha=ha, va='center',
                fontsize=6.5, color=INK_2, zorder=5)


def layout_handle(layout, mode):
    color = LAYOUT_COLOR[layout]
    return Line2D([], [], linestyle='none', marker=MODE_MARKER[mode], markersize=5.5, markerfacecolor=color,
                  markeredgecolor=color, markeredgewidth=0.5, label=LAYOUT_NAME[layout] + MODE_SUFFIX[mode])


def rope_ms_handle(length):
    if length == 1:
        return Line2D([], [], linestyle='none', marker='o', markersize=5.5, markerfacecolor=ROPE_RED,
                      markeredgecolor=SURFACE, markeredgewidth=0.5, label='ropebwt3 -l1')
    color, marker = ROPE_LONG_STYLE[length]
    return Line2D([], [], linestyle='none', marker=marker, markersize=5.5, markerfacecolor=color,
                  markeredgecolor=SURFACE, markeredgewidth=0.5, label=f'ropebwt3 -l{length}*')


def rope_shift_handle(shift):
    return Line2D([], [], linestyle='-', linewidth=0.8, color=ROPE_RED, marker=ROPE_SHIFT_MARKER[shift],
                  markersize=6.5 if shift == 8 else 5.5, markerfacecolor=ROPE_RED, markeredgecolor=SURFACE,
                  markeredgewidth=0.5, label=f'ropebwt3 -s{shift}')


def others_handle(name, marker):
    return Line2D([], [], linestyle='none', marker=marker, markersize=5.5, markerfacecolor=SURFACE,
                  markeredgecolor=OTHERS_COLOR, markeredgewidth=1.2, label=name)


def mode_key_handles(modes):
    """Gray marker key so the mode encoding reads independently of layout color."""
    out = []
    for mode in modes:
        if mode in REPLAY_TICK:
            size, width = REPLAY_TICK[mode]
            out.append(Line2D([], [], linestyle='none', marker='_', markersize=size, markeredgewidth=width,
                              color=INK_2, label=mode))
        else:
            out.append(Line2D([], [], linestyle='none', marker=MODE_MARKER[mode], markersize=5.5,
                              markerfacecolor=INK_2, markeredgecolor=INK_2, markeredgewidth=0.5, label=mode))
    return out


BLANK = Line2D([], [], linestyle='none', label=' ')


def pad_group(handles, rows=2):
    """Pads a legend group to whole columns so groups never share a column."""
    return handles + [BLANK] * (-len(handles) % rows)


def save(fig, out_dir, stem):
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(out_dir / f'{stem}.{ext}', dpi=300, facecolor=SURFACE, bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)


def legend_below(fig, handles, ncol, y=0.0):
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, y), ncol=ncol, frameon=False,
               handletextpad=0.3, columnspacing=1.1)


# MS figures

MS_YLABEL = 'Mbases/s (query time)'
MEM_XLABEL = 'Peak resident memory (GiB)'


def draw_ms_panel(ax, cohort, summary, rec, figure, replay=False, ftab=False, others=False, label_points=True,
                  layouts=LAYOUTS):
    cs = cells(summary)
    for layout in layouts:
        ordinary = cs.get(f'{layout}-ordinary')
        if not ok(ordinary):
            continue
        color = LAYOUT_COLOR[layout]
        x0, y0 = gib(ordinary), ordinary['mbases_per_second']
        point(ax, x0, y0, color, MODE_MARKER['ordinary'], z=5)
        rec.add(figure, cohort, LAYOUT_NAME[layout], x0, y0, 'Mbases/s', ordinary)
        if replay:
            got = [(m, cs.get(f'{layout}-{m}')) for m in REPLAY_MODES]
            got = [(m, c) for m, c in got if ok(c)]
            if got:
                ys = [y0] + [c['mbases_per_second'] for _, c in got]
                ax.plot([x0, x0], [min(ys), max(ys)], color=color, linewidth=1.0, solid_capstyle='butt',
                        zorder=2)
                for mode, c in got:
                    y = c['mbases_per_second']
                    replay_tick(ax, x0, y, color, mode)
                    rec.add(figure, cohort, LAYOUT_NAME[layout] + MODE_SUFFIX[mode], x0, y, 'Mbases/s', c,
                            x_note='drawn as a tick at the ordinary peak memory')
        if ftab:
            c = cs.get(f'{layout}-ftab')
            if ok(c):
                x, y = gib(c), c['mbases_per_second']
                ax.annotate('', xy=(x, y), xytext=(x0, y0), zorder=2,
                            arrowprops=dict(arrowstyle='-', color=color, linewidth=0.8, alpha=0.6,
                                            shrinkA=3, shrinkB=3))
                point(ax, x, y, color, MODE_MARKER['ftab'], z=4)
                rec.add(figure, cohort, LAYOUT_NAME[layout] + '+ftab', x, y, 'Mbases/s', c)
        if label_points:
            direct_label(ax, x0, y0, LAYOUT_NAME[layout].replace('ms2-', ''), dx=-5, ha='right')
    for cid, length in ROPE_MS:
        c = cs.get(cid)
        if not ok(c):
            continue
        x, y = gib(c), c['mbases_per_second']
        if length == 1:
            point(ax, x, y, ROPE_RED, 'o', z=4)
            if label_points:
                direct_label(ax, x, y, '-l1', dx=6)
            label = 'ropebwt3 -l1'
        else:
            color, marker = ROPE_LONG_STYLE[length]
            point(ax, x, y, color, marker, z=3)
            label = f'ropebwt3 -l{length}*'
        rec.add(figure, cohort, label, x, y, 'Mbases/s', c)
    if others:
        for cid, name, marker, dx, dy, ha in OTHERS:
            c = cs.get(cid)
            if not ok(c):
                continue
            x, y = gib(c), c['mbases_per_second']
            point(ax, x, y, OTHERS_COLOR, marker, hollow=True, z=4)
            rec.add(figure, cohort, name, x, y, 'Mbases/s', c, x_note='one-read-wall fallback timing')
            if label_points:
                direct_label(ax, x, y, name, dx=dx, dy=dy, ha=ha)


def ms_limits(ax, summaries, ids):
    pts = [cells(s).get(i) for s in summaries for i in ids]
    pts = [c for c in pts if ok(c)]
    ys = [c['mbases_per_second'] for c in pts]
    xs = [gib(c) for c in pts]
    ax.set_ylim(min(ys) / 1.6, max(ys) * 1.8)
    ax.set_xlim(0, max(xs) * 1.08)
    ax.yaxis.set_major_locator(FixedLocator([0.01, 0.1, 1, 10, 100]))
    ax.yaxis.set_major_formatter(decade_formatter())


def ms_grid_figure(ms, stem, out_dir, rec, replay):
    present = [c for c in COHORTS if c in ms]
    fig, axes = plt.subplots(2, 2, figsize=(FIG_W, 5.0), sharex=True, sharey=True, facecolor=SURFACE)
    for ax, cohort in zip(axes.flat, present):
        draw_ms_panel(ax, cohort, ms[cohort], rec, stem, replay=replay, layouts=MS_GRID_LAYOUTS)
        ax.set_title(TITLES[cohort], loc='left')
    for ax in axes.flat:
        style_axes(ax, MS_YLABEL, MEM_XLABEL)
        ms_limits(ax, [ms[c] for c in present], HUMAN_MS_IDS)
        ax.label_outer()
    # The legend fills column by column, two rows deep; padding keeps each group in its own columns.
    handles = pad_group([layout_handle(l, 'ordinary') for l in MS_GRID_LAYOUTS])
    if replay:
        handles += mode_key_handles(['ordinary', 'oracle'])
    handles += pad_group([rope_ms_handle(length) for _, length in ROPE_MS])
    fig.tight_layout(h_pad=1.2, w_pad=1.0)
    legend_below(fig, handles, ncol=len(handles) // 2)
    save(fig, out_dir, stem)


def ms_ftab_figure(ms, out_dir, rec):
    stem = 'fig-ms-ftab'
    summary = ms['human_short']
    fig, ax = plt.subplots(figsize=(3.45, 3.0), facecolor=SURFACE)
    draw_ms_panel(ax, 'human_short', summary, rec, stem, ftab=True, layouts=MS_GRID_LAYOUTS)
    ax.set_title(TITLES['human_short'], loc='left')
    style_axes(ax, MS_YLABEL, MEM_XLABEL)
    ms_limits(ax, [summary], HUMAN_MS_IDS)
    handles = [layout_handle(l, 'ordinary') for l in MS_GRID_LAYOUTS if f'{l}-ordinary' in cells(summary)]
    handles += mode_key_handles(['ordinary', 'ftab'])
    handles += [rope_ms_handle(length) for _, length in ROPE_MS]
    fig.tight_layout()
    legend_below(fig, handles, ncol=2)
    save(fig, out_dir, stem)


def chr22_figure(chr22, h, out_dir, rec):
    stem = f'fig-ms22-{h}'
    present = [s for s in CHR22_SETS if (h, s) in chr22]
    if not present:
        return False
    sums = [chr22[(h, s)] for s in present]
    layouts = ['ms', 'psi', 'phi', 'dual'] if h in CHR22_PHI_REFS else LAYOUTS
    fig, axes = plt.subplots(1, len(present), figsize=(FIG_W, 2.75), sharex=True, sharey=True,
                             facecolor=SURFACE, squeeze=False)
    for ax, s in zip(axes.flat, present):
        draw_ms_panel(ax, s, chr22[(h, s)], rec, stem, others=True, layouts=layouts)
        ax.set_title(TITLES[s], loc='left')
    for ax in axes.flat:
        style_axes(ax, MS_YLABEL, MEM_XLABEL)
        ms_limits(ax, sums, CHR22_MS_IDS)
        ax.label_outer()
    handles = pad_group([layout_handle(l, 'ordinary') for l in layouts])
    handles += pad_group([rope_ms_handle(length) for _, length in ROPE_MS])
    handles += pad_group([others_handle(name, marker) for _, name, marker, *_ in OTHERS])
    fig.tight_layout(w_pad=1.0)
    legend_below(fig, handles, ncol=len(handles) // 2)
    save(fig, out_dir, stem)
    return True


# Locate figures

LOC_YLABEL = 'Reported locations/s (query time)'
LOC_LAYOUTS = ['psi', 'phi', 'dual']
# Psi, Phi and Dual locate from one index file, so they share a peak; a small
# horizontal offset keeps each visible. The CSV keeps the measured value.
DODGE = {'psi': -8.0, 'phi': 0.0, 'dual': 8.0}
# With ftab the six ms2 points sit within 4 GiB of each other; spread them
# wider (psi, psi+ftab, phi, phi+ftab, dual, dual+ftab from left to right).
FTAB_DODGE = {('psi', False): -25.0, ('psi', True): -15.0, ('phi', False): -5.0, ('phi', True): 5.0,
              ('dual', False): 15.0, ('dual', True): 25.0}


def draw_locate_panel(ax, cohort, summary, cap, rec, figure, ftab=False, label_points=True, replay=False):
    """With replay, each ms2 layout also gets its oracle replay as a tick on a vertical line above the
ordinary point (oracle-scan is not drawn)."""
    cs = cells(summary)
    for layout in LOC_LAYOUTS:
        color = LAYOUT_COLOR[layout]
        c = cs.get(f'cap{cap}-{layout}-ordinary')
        if not ok(c):
            continue
        base = gib(c)
        x0 = base + (FTAB_DODGE[(layout, False)] if ftab else DODGE[layout])
        y0 = c['locations_per_second']
        point(ax, x0, y0, color, MODE_MARKER['ordinary'], z=4)
        rec.add(figure, cohort, LAYOUT_NAME[layout], x0, y0, 'locations/s', c,
                x_note='offset horizontally for visibility')
        o = cs.get(f'cap{cap}-{layout}-oracle') if replay else None
        if ok(o):
            yo = o['locations_per_second']
            ax.plot([x0, x0], [min(y0, yo), max(y0, yo)], color=color, linewidth=1.0, solid_capstyle='butt',
                    zorder=2)
            replay_tick(ax, x0, yo, color, 'oracle')
            rec.add(figure, cohort, LAYOUT_NAME[layout] + MODE_SUFFIX['oracle'], x0, yo, 'locations/s', o,
                    x_note='drawn as a tick above the ordinary point')
        if ftab:
            f = cs.get(f'cap{cap}-{layout}-ftab')
            if ok(f):
                x, y = base + FTAB_DODGE[(layout, True)], f['locations_per_second']
                ax.annotate('', xy=(x, y), xytext=(x0, y0), zorder=2,
                            arrowprops=dict(arrowstyle='-', color=color, linewidth=0.8, alpha=0.6,
                                            shrinkA=3, shrinkB=3))
                point(ax, x, y, color, MODE_MARKER['ftab'], z=4)
                rec.add(figure, cohort, LAYOUT_NAME[layout] + '+ftab', x, y, 'locations/s', f,
                        x_note='offset horizontally for visibility')
                direct_label(ax, x, y, f'×{y / y0:.2f}', dx=0, dy=8 if layout == 'psi' else -8, ha='center')
    if label_points:
        pts = [cs.get(f'cap{cap}-{l}-ordinary') for l in LOC_LAYOUTS]
        if all(ok(p) for p in pts):
            hi = max(pts, key=lambda p: p['locations_per_second'])
            direct_label(ax, gib(hi) - (24 if ftab else 9), hi['locations_per_second'], 'psi, phi, dual', dx=-3,
                         ha='right')
    xs, ys = [], []
    for shift in ROPE_SHIFTS:
        c = cs.get(f'cap{cap}-rope-s{shift}')
        if not ok(c):
            continue
        x, y = gib(c), c['locations_per_second']
        xs.append(x)
        ys.append(y)
        point(ax, x, y, ROPE_RED, ROPE_SHIFT_MARKER[shift], size=MARKER_SIZE * (1.6 if shift == 8 else 1.1), z=4)
        rec.add(figure, cohort, f'ropebwt3 -s{shift}', x, y, 'locations/s', c)
        if label_points:
            direct_label(ax, x, y, f'-s{shift}', dx=0, dy=-9, ha='center')
    if len(xs) > 1:
        ax.plot(xs, ys, color=ROPE_RED, linewidth=0.8, alpha=0.6, zorder=2)


def locate_limits(ax, summaries):
    ys = [c['locations_per_second'] for s in summaries for c in s['cells'] if ok(c)]
    xs = [gib(c) for s in summaries for c in s['cells'] if ok(c)]
    ax.set_ylim(min(ys) / 2.2, max(ys) * 2.0)
    ax.set_xlim(0, max(xs) * 1.08)
    ax.yaxis.set_major_locator(FixedLocator([1e3, 1e4, 1e5, 1e6, 1e7]))
    ax.yaxis.set_major_formatter(decade_formatter())


def locate_handles(ftab=False):
    handles = [layout_handle(l, 'ordinary') for l in LOC_LAYOUTS]
    if ftab:
        handles += mode_key_handles(['ordinary', 'ftab'])
    handles += [rope_shift_handle(s) for s in ROPE_SHIFTS]
    return handles


def locate_grid_figure(loc, cap, out_dir, rec):
    stem = f'fig-l{cap}'
    present = [c for c in COHORTS if c in loc]
    fig, axes = plt.subplots(2, 2, figsize=(FIG_W, 5.0), sharex=True, sharey=True, facecolor=SURFACE)
    for ax, cohort in zip(axes.flat, present):
        draw_locate_panel(ax, cohort, loc[cohort], cap, rec, stem)
        ax.set_title(TITLES[cohort], loc='left')
    for ax in axes.flat:
        style_axes(ax, LOC_YLABEL, MEM_XLABEL)
        locate_limits(ax, [loc[c] for c in present])
        ax.label_outer()
    fig.tight_layout(h_pad=1.2, w_pad=1.0)
    legend_below(fig, locate_handles(), ncol=5)
    save(fig, out_dir, stem)


def locate_oracle_figure(loc, cap, out_dir, rec):
    stem = LORACLE_CAPS[cap]
    present = [c for c in COHORTS if c in loc]
    fig, axes = plt.subplots(2, 2, figsize=(FIG_W, 5.0), sharex=True, sharey=True, facecolor=SURFACE)
    for ax, cohort in zip(axes.flat, present):
        draw_locate_panel(ax, cohort, loc[cohort], cap, rec, stem, replay=True)
        ax.set_title(TITLES[cohort], loc='left')
    for ax in axes.flat:
        style_axes(ax, LOC_YLABEL, MEM_XLABEL)
        locate_limits(ax, [loc[c] for c in present])
        ax.label_outer()
    fig.tight_layout(h_pad=1.2, w_pad=1.0)
    legend_below(fig, locate_handles() + mode_key_handles(['oracle']), ncol=6)
    save(fig, out_dir, stem)


def oracle_speedups(loc, cap):
    """{(cohort, layout): oracle over ordinary locations/s} at one cap."""
    out = {}
    for cohort in COHORTS:
        cs = cells(loc[cohort]) if cohort in loc else {}
        for layout in LOC_LAYOUTS:
            a, b = cs.get(f'cap{cap}-{layout}-ordinary'), cs.get(f'cap{cap}-{layout}-oracle')
            if ok(a) and ok(b):
                out[(cohort, layout)] = b['locations_per_second'] / a['locations_per_second']
    return out


def locate_ftab_figure(loc, out_dir, rec):
    stem = 'fig-lftab'
    summary = loc['human_short']
    fig, axes = plt.subplots(1, 3, figsize=(FIG_W, 2.6), sharex=True, sharey=True, facecolor=SURFACE)
    for ax, cap in zip(axes, CAPS):
        draw_locate_panel(ax, 'human_short', summary, cap, rec, stem, ftab=True)
        ax.set_title(f'N = {cap}', loc='left')
        style_axes(ax, LOC_YLABEL, MEM_XLABEL)
        locate_limits(ax, [summary])
        ax.label_outer()
    fig.tight_layout(w_pad=1.0)
    legend_below(fig, locate_handles(ftab=True), ncol=7)
    save(fig, out_dir, stem)


# Ops figure: time and operation counts per base on held-out genomes

OPS_LAYOUTS = ['ms', 'psi', 'dual', 'phi']
OPS_NAME = {'ms': 'lcp', 'psi': 'psi', 'dual': 'dual', 'phi': 'phi'}
TIME_YLABEL = 'ns per base (query time)'
OPS_YLABEL = 'Operations per base'
BAR_W = 0.64
# Light tint of a layout color for the LF segment of an ops bar.
LF_TINT = 0.55


def tint(hex_color, amount):
    r, g, b = (int(hex_color[i:i + 2], 16) for i in (1, 3, 5))
    mix = [round(c + (255 - c) * amount) for c in (r, g, b)]
    return '#' + ''.join(f'{c:02x}' for c in mix)


def ops_breakdown(cell):
    """Returns (LF steps, other steps) per base; other is psi plus phi plus walk rows."""
    o = cell['ops_per_base']
    return o['lf'], o['psi'] + o['phi'] + o['walk']


def bar_label(ax, x, y, text):
    ax.annotate(text, (x, y), xytext=(0, 1.5), textcoords='offset points', ha='center', va='bottom',
                fontsize=6.3, color=INK_2, zorder=6)


def ops_axes(ax, ylabel, show_ylabel=True):
    ax.grid(True, axis='y', color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.set_xticks(range(len(OPS_LAYOUTS)), [OPS_NAME[l] for l in OPS_LAYOUTS])
    ax.tick_params(axis='x', length=0)
    ax.set_xlim(-0.6, len(OPS_LAYOUTS) - 0.4)
    if show_ylabel:
        ax.set_ylabel(ylabel)


def draw_time_panel(ax, ref, summary, rec, stem='fig-ops'):
    cs = ops_cells(summary)
    top = 0.0
    for i, layout in enumerate(OPS_LAYOUTS):
        c = cs.get(f'{layout}-ordinary')
        if c is not None and c.get('status') == 'capped':
            # A capped cell (Fig Ops2 ms2-phi on human reads) has no time; its label sits on the axis.
            bar_label(ax, i, 0, 'capped')
        if not ok(c):
            continue
        y = c['ns_per_base']
        ax.bar(i, y, width=BAR_W, color=LAYOUT_COLOR_OPS[layout], edgecolor=SURFACE, linewidth=0.8, zorder=3)
        bar_label(ax, i, y, f'{y:.0f}')
        rec.add(stem, ref, f'ms2-{OPS_NAME[layout]}', gib(c), y, 'ns/base', c)
        top = max(top, y)
        r = cs.get(f'{layout}-oracle')
        if ok(r):
            yr = r['ns_per_base']
            half = BAR_W * 0.66
            ax.plot([i - half, i + half], [yr, yr], color=INK, linewidth=1.3, solid_capstyle='butt', zorder=4)
            rec.add(stem, ref, f'ms2-{OPS_NAME[layout]}{MODE_SUFFIX["oracle"]}', gib(r), yr, 'ns/base', r,
                    x_note='drawn as a tick across the ordinary bar')
    ax.set_ylim(0, (top or 1.0) * 1.16)


def draw_ops_panel(ax, ref, summary, rec, stem='fig-ops'):
    """Stacked operation bars; returns the tallest total."""
    cs = ops_cells(summary)
    top = 0.0
    for i, layout in enumerate(OPS_LAYOUTS):
        c = cs.get(f'{layout}-ordinary')
        if not ok(c):
            continue
        lf, other = ops_breakdown(c)
        color = LAYOUT_COLOR_OPS[layout]
        # The LF segment is what the oracle performs: replay takes one LF step per base and nothing else.
        ax.bar(i, lf, width=BAR_W, color=tint(color, LF_TINT), edgecolor=SURFACE, linewidth=0.8, zorder=3)
        if other > 0:
            ax.bar(i, other, bottom=lf, width=BAR_W, color=color, edgecolor=SURFACE, linewidth=0.8, zorder=3)
        bar_label(ax, i, lf + other, f'{lf + other:.2f}')
        rec.add(stem, ref, f'ms2-{OPS_NAME[layout]} ops', gib(c), lf + other, 'ops/base', c,
                x_note=f'LF {lf:.4f} + psi, phi, walk {other:.4f}')
        top = max(top, lf + other)
    return top


def ops_figure(ops, out_dir, rec, stem='fig-ops'):
    present = [r for r in OPS_REFS if r in ops]
    if not present:
        return False
    fig = plt.figure(figsize=(FIG_W, 2.05), facecolor=SURFACE)
    gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 0.16, 1, 1], wspace=0.42)
    time_axes = [fig.add_subplot(gs[0, i]) for i in (0, 1)]
    ops_axes_ = [fig.add_subplot(gs[0, 3])]
    ops_axes_.append(fig.add_subplot(gs[0, 4], sharey=ops_axes_[0]))
    top = 0.0
    for k, ref in enumerate(OPS_REFS):
        if ref not in ops:
            continue
        draw_time_panel(time_axes[k], ref, ops[ref], rec, stem)
        top = max(top, draw_ops_panel(ops_axes_[k], ref, ops[ref], rec, stem))
        ops_axes(time_axes[k], TIME_YLABEL, show_ylabel=(k == 0))
        ops_axes(ops_axes_[k], OPS_YLABEL, show_ylabel=(k == 0))
        time_axes[k].set_title(OPS_TITLE[ref], loc='left')
        ops_axes_[k].set_title(OPS_TITLE[ref], loc='left')
    # Fig Ops totals stay under 4 operations per base; read queries can reposition far more often.
    step = 1 if top <= 4 else next((s for s in (1, 2, 5, 10, 20, 50, 100) if top / s <= 5), 200)
    ops_max = 4 if top <= 4 else step * -(-top // step)
    for ax in ops_axes_:
        ax.set_ylim(0, ops_max)
        ax.yaxis.set_major_locator(FixedLocator(list(range(0, int(ops_max) + 1, step))))
    plt.setp(ops_axes_[1].get_yticklabels(), visible=False)
    handles = [Patch(facecolor=LAYOUT_COLOR_OPS[l], edgecolor='none', label=OPS_NAME[l]) for l in OPS_LAYOUTS]
    handles += [Line2D([], [], color=INK, linewidth=1.3, label='oracle time'),
                Patch(facecolor=tint(INK_2, LF_TINT), edgecolor='none', label='LF steps (oracle\'s ops)'),
                Patch(facecolor=INK_2, edgecolor='none', label='added psi, phi, walk steps')]
    fig.subplots_adjust(left=0.07, right=0.995, top=0.86, bottom=0.2)
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.01), ncol=len(handles),
               frameon=False, handlelength=1.4, handletextpad=0.35, columnspacing=0.9)
    save(fig, out_dir, stem)
    return True


# Captions

def fmt_int(n):
    return f'{n:,}'


def fmt_mb(bases):
    return f'{bases / 1e6:.2f} Mbases' if bases >= 1e6 else f'{bases / 1e3:.0f} kbases'


def query_sets(summaries, cohorts):
    parts = []
    for cohort in cohorts:
        if cohort in summaries:
            reads, bases = workload(summaries[cohort])
            parts.append(f'{TITLES[cohort]}, {fmt_int(reads)} reads, {fmt_mb(bases)}')
    return '; '.join(parts)


def common_tail(summaries, reference=HUMAN472):
    revs = ', '.join(revisions(summaries))
    return (f'Reference: {reference}. Single thread, NUMA-bound (numactl CPU and memory '
            f'binding to one node). Speed is query time only; index load time is excluded and reported separately '
            f'(Load times table below). Each value is the median of two timed runs after one warmup run. Orbit '
            f'commit {revs}. Peak resident memory is the process peak RSS.')


def load_table(rows, figure):
    sub = [r for r in rows if r['figure'] == figure]
    seen = {}
    for r in sub:
        key = (r['read_set'], r['series'])
        if key not in seen:
            seen[key] = r
    lines = ['| Read set | Series | Median load (s) | Median query (s) |', '| --- | --- | ---: | ---: |']
    for (cohort, series), r in seen.items():
        load = r['median_load_s']
        lines.append(f'| {TITLES[cohort]} | {series} | {load:.2f} | {r["median_query_s"]:.2f} |'
                     if load is not None else f'| {TITLES[cohort]} | {series} | n/a | n/a |')
    return '\n'.join(lines)


def missing_cells(ms, chr22, loc):
    out = []
    for cohort in COHORTS:
        if cohort not in ms:
            out.append(f'MS summary for {TITLES[cohort]} not found.')
            continue
        cs = cells(ms[cohort])
        for cid in HUMAN_MS_IDS:
            # ftab cells exist only for short reads, the one set Fig MS+ftab draws.
            if cid.endswith('-ftab') and cohort != 'human_short':
                continue
            c = cs.get(cid)
            if not ok(c):
                out.append(f'MS {TITLES[cohort]}: {cid} {"not run" if c is None else c["status"]}.')
    for h in CHR22_REFS:
        for s in CHR22_SETS:
            if (h, s) not in chr22:
                continue
            cs = cells(chr22[(h, s)])
            for cid in CHR22_MS_IDS:
                c = cs.get(cid)
                if not ok(c):
                    out.append(f'MS chr22-h{h} {TITLES[s]}: {cid} {"not run" if c is None else c["status"]}.')
    for cohort in COHORTS:
        if cohort not in loc:
            out.append(f'Locate summary for {TITLES[cohort]} not found.')
            continue
        bad = [c['id'] for c in loc[cohort]['cells'] if not ok(c)]
        if bad:
            out.append(f'Locate {TITLES[cohort]}: not checked: {", ".join(bad)}.')
    return out


def chr22_caption(chr22, h):
    present = [s for s in CHR22_SETS if (h, s) in chr22]
    missing = [s for s in CHR22_SETS if (h, s) not in chr22]
    sets = []
    for s in CHR22_SETS:
        if (h, s) in chr22:
            reads, bases = workload(chr22[(h, s)])
        else:
            reads, bases = CHR22_SIM_SIZE[s]
        sets.append(f'{CHR22_SIM_TEXT[s]}; {fmt_int(reads)} reads, {fmt_mb(bases)}.')
    rope_k = {chr22[(h, s)]['config'].get('rope', {}).get('batch_size') for s in present}
    k_note = ' ropebwt3 runs with -K 1m.' if rope_k == {'1m'} else ''
    body = (f'Matching statistics (MS) only on {CHR22_REF_TEXT[h].split(" (")[0]}, simulated reads: query speed '
            'against peak resident memory for ms2 layouts (ms2-lcp, ms2-psi, ms2-dual, ordinary queries), '
            'ropebwt3 MS-only modes, and MONI, PHONI and AugPHONI. Reads are simulated from chr22 of the held-out '
            'haplotype NA19240 paternal (HPRC sample 200245_NA19240.pat; chr22 contigs placed by its RagTag '
            'scaffold, 51.09 Mbases of ACGT), which is absent from both the chr22-h10 and the chr22-h100 '
            'references (neither the sample nor its donor is in either panel). Each set targets 60 Mbases of '
            'whole ACGT-only reads. ' + ' '.join(sets) +
            ' Manifests: results/campaign-20261012/inputs/chr22_sim/ (short-60mb-v1, pacbio-60mb-v1, '
            f'ont-60mb-v1, haplotype-NA19240-pat-chr22-v1). ropebwt3 -l1 (red) reports every MS value. '
            f'{ASTERISK_NOTE}{k_note} {OTHERS_NOTE}')
    if missing:
        body += (' Pending: the ' + ', '.join(TITLES[s] for s in missing) +
                 ' panel is not drawn because its run had not finished when this draft was made.')
    return body, common_tail([chr22[(h, s)] for s in present], reference=CHR22_REF_TEXT[h])


def ops_caption(ops):
    sizes = {r: workload(ops[r]) for r in ops}
    rep = {r: ops_cells(ops[r])['ms-ordinary']['ops_per_base']['repositions'] for r in ops}
    m_reads, m_bases = sizes.get('mtb152', (16, 70570430))
    h_reads, h_bases = sizes.get('human472', (81, 3038047138))
    return ('Time and operations per base for held-out genomes: ms2-lcp, ms2-psi, ms2-dual and ms2-phi '
            '(MS only, ordinary queries). Left pair: nanoseconds per query base, one bar per layout; a black '
            'tick across the bar marks the same layout\'s oracle replay. Right '
            'pair: operations per query base, stacked as LF steps (light segment) and the steps the layout adds '
            'on top (full color): psi steps for ms2-psi, phi steps for ms2-phi, both for ms2-dual, and walk '
            'rows (LCP rows scanned while repositioning) for ms2-lcp. The oracle takes exactly one LF step per '
            'base and nothing else, so the light segment is also the oracle\'s operation count; the full-color '
            'segment is the work the oracle avoids. Queries are genomes absent from the index. mtb152 is the '
            'TeraTools paper\'s 95 percent split queried with the 8 held-out genomes, each followed by its '
            f'reverse complement ({fmt_int(m_reads)} records, {m_bases / 1e6:.1f} Mbases). HPRC2 full is '
            'Human472 queried with APR008 haplotype 1 from the Arab Pangenome Reference (GenBank '
            'GCA_050494745.1, OpenHGL 210019_APR008.hap1), absent from Human472 '
            f'({fmt_int(h_reads)} ACGT-only records, {h_bases / 1e9:.2f} Gbases). Records are split at every '
            'run of characters other than A, C, G and T. Repositions (mismatches that force the query to '
            f'restart its match) are rare: {rep.get("mtb152", 0):.5f} per base on mtb152 and '
            f'{rep.get("human472", 0):.5f} on HPRC2 full, which is why replay saves so little for ms2-lcp, '
            'whose only extra work is the walk at a reposition.')


def ops2_caption(ops2, reads=OPS2_DEFAULT_READS):
    sizes = {r: workload(ops2[r]) for r in ops2}
    rep = {r: (ops_cells(ops2[r]).get('ms-ordinary') or {}).get('ops_per_base', {}).get('repositions') for r in ops2}
    parts = [f'{OPS_TITLE[r]} {fmt_int(sizes[r][0])} reads, {sizes[r][1] / 1e6:.1f} Mbases' for r in OPS_REFS if r in sizes]
    reps = [f'{rep[r]:.4f} on {OPS_TITLE[r]}' for r in OPS_REFS if rep.get(r) is not None]
    capped = [f'ms2-{OPS_NAME[c["layout"]]} {c["mode"]} on {OPS_TITLE[r]} (cap {c.get("cap_s")} s)'
              for r in OPS_REFS if r in ops2 for c in ops2[r]['cells']
              if c.get('interleave') == 32 and c.get('status') == 'capped' and c['mode'] != 'oracle-scan']
    text = ('As Fig Ops, with ' + OPS2_SIM_TEXT[reads] + ' mtb152 reads come from the forward strand of the 8 '
            'held-out genomes and HPRC2 full reads from APR008 haplotype 1 '
            f'({"; ".join(parts)}).')
    if reps:
        text += f' Repositions per base: {", ".join(reps)}.'
    if capped:
        text += (' Capped (a per-command time limit was reached, so no bar is drawn): ' + ', '.join(capped) + '.')
    return text


def ops_tail(ops):
    summaries = [ops[r] for r in OPS_REFS if r in ops]
    revs = ', '.join(revisions(summaries))
    return (f'References: mtb152 and {HUMAN472}. Single thread, interleave 32 (queries in flight at once), '
            'NUMA-bound (numactl CPU and memory binding to one node). Time is query time only; index load time '
            'is excluded and reported separately (Load times table below). Each time is the median of two '
            'timed runs after one warmup run. Operation counts come from a counting build of the same commit; '
            f'Orbit commit {revs}.')


def write_captions(ms, chr22, loc, ops, rows, chr22_drawn, out_dir, ops2=None, ops2_reads=OPS2_DEFAULT_READS):
    ms_all = [ms[c] for c in COHORTS if c in ms]
    loc_all = [loc[c] for c in COHORTS if c in loc]
    ms_sets = query_sets(ms, COHORTS)
    loc_sets = query_sets(loc, COHORTS)
    ms_tail = common_tail(ms_all)
    loc_tail = common_tail(loc_all)
    loc_k = {s['config'].get('rope', {}).get('batch_size') for s in loc_all}
    loc_scope = ('ms2-min has no locate support, so the ms2 SMEM reporting set is ms2-psi and ms2-dual '
                 '(ordinary). Both locate from the same index file (the Dual index with sequence ids), so they '
                 'share one peak memory. ropebwt3 runs -l1 -pN with sampled suffix arrays at shifts 11, 8 and 6 '
                 '(-s11, -s8, -s6; denser sampling uses more memory)' +
                 (', with -K 1m' if loc_k == {'1m'} else '') + '. Every output was checked: ms2 and ropebwt3 '
                 'agree with zero differences at each cap.')
    out = ['# Draft figure captions', '',
           'Generated by `bench/campaign-20261012/paper_figures.py`. Rerun the script to '
           'refresh the figures, captions and `plotted-points.csv` after the frozen runs.', '']

    def fig(stem, title, body, tail):
        out.extend([f'## {title} (`{stem}`)', '', f'![{title}]({stem}.png)', '', body + ' ' + tail, '',
                    '<details><summary>Load times</summary>', '', load_table(rows, stem), '', '</details>', ''])

    fig('fig-ms', 'Fig MS',
        'Matching statistics (MS) only: query speed against peak resident memory for ms2 layouts '
        '(ms2-lcp, ms2-psi, ms2-dual, ordinary queries) and ropebwt3 MS-only modes. Query sets: '
        f'{ms_sets}; all reads are whole ACGT-only records. ropebwt3 -l1 (red) reports every MS value. '
        f'{ASTERISK_NOTE}', ms_tail)
    fig('fig-ms-oracle', 'Fig MS+Oracle',
        'As Fig MS, adding oracle replay for each ms2 layout. A vertical line in the layout\'s color rises from '
        f'the ordinary point (circle) to a horizontal tick for the oracle. {ORACLE_NOTE} '
        f'Query sets: {ms_sets}. {ASTERISK_NOTE}', ms_tail)
    fig('fig-ms-ftab', 'Fig MS+ftab',
        'Short reads only: MS speed against peak memory for ms2 layouts without (circle) and with (diamond) an '
        'ftab (k-mer lookup table that skips the first steps of backward search), plus ropebwt3 MS-only modes. '
        'Lines join each layout with its ftab variant; the ftab point is at its own measured peak. '
        f'Query set: {query_sets(ms, ["human_short"])}. {ASTERISK_NOTE}',
        common_tail([ms['human_short']]))
    for h in CHR22_REFS:
        if h in chr22_drawn:
            body, tail = chr22_caption(chr22, h)
            fig(f'fig-ms22-{h}', f'Fig MS22-{h}', body, tail)
        else:
            out.extend([f'## Fig MS22-{h} (`fig-ms22-{h}`)', '', 'Not drawn: no chr22 summary found.', ''])
    for cap in CAPS:
        fig(f'fig-l{cap}', f'Fig L{cap}',
            f'SMEM offset reporting (locate) with at most N = {cap} reported positions per SMEM: reported locations '
            f'per second of query time against peak memory. {loc_scope} ms2-psi and ms2-dual are offset '
            f'horizontally by 5 GiB to either side of their shared peak for visibility. Query sets (smaller than '
            f'the MS sets): {loc_sets}.', loc_tail)
    fig('fig-lftab', 'Fig LFtab',
        'Short reads only, SMEM offset reporting at N = 1, 10 and 100: ms2-psi and ms2-dual without (circle) and '
        'with (diamond) ftab, plus ropebwt3 -s11, -s8 and -s6. Each ftab point is labeled with its speed ratio '
        'over the same layout without ftab. The four ms2 points are spread horizontally for visibility (left '
        'to right: psi, psi+ftab, dual, dual+ftab); measured peaks are about 293 GiB without ftab and 296 GiB '
        'with it. ' + loc_scope +
        f' Query set: {query_sets(loc, ["human_short"])}.', common_tail([loc['human_short']]))
    for cap, stem in sorted(LORACLE_CAPS.items(), key=lambda item: -item[0]):
        if loc:
            fig(stem, f'Fig LOracle{"" if cap == 100 else cap}',
                f'SMEM offset reporting (locate) with at most N = {cap} reported positions per SMEM, as Fig L{cap}, '
                'adding oracle replay for ms2-psi and ms2-dual. A vertical line in the layout\'s color rises from '
                'the ordinary point (circle) to a horizontal tick for the oracle, which replays a recorded run with '
                'every repositioning decision known in advance; its output is identical to the ordinary run\'s. '
                'The oracle replaces only the matching-statistics phase; counting, listing and formatting the '
                'reported locations run unchanged. ' + loc_scope + ' ms2-psi and ms2-dual are offset horizontally '
                'by 5 GiB to either side of their shared peak for visibility. Oracle over ordinary locations per '
                'second: ' + '; '.join(f'{TITLES[c]} {LAYOUT_NAME[l].replace("ms2-", "")} ×{r:.2f}'
                                         for (c, l), r in oracle_speedups(loc, cap).items()) +
                f'. Query sets (smaller than the MS sets): {loc_sets}.', loc_tail)
    if ops:
        fig('fig-ops', 'Fig Ops', ops_caption(ops), ops_tail(ops))
    if ops2:
        fig('fig-ops2', 'Fig Ops2', ops2_caption(ops2, ops2_reads), ops_tail(ops2))

    still = [
        '- Freeze reruns with ropebwt3 -K 1m for the Human472 MS cells, which used the default -K 100m. At -l1 '
        'that buffers SMEMs, visible on the skin metagenome (15.9 GiB at -l1 against 14.5 GiB for other '
        'lengths). The locate figures already use the -K 1m oracle-locate runs.',
    ]
    for h in CHR22_REFS:
        gone = [TITLES[s] for s in CHR22_SETS if (h, s) not in chr22]
        if gone:
            still.append(f'- Fig MS22-{h}: {", ".join(gone)} panel pending (summary not yet on Rockfish at '
                         f'{CHR22_REMOTE.format(h=h, s="simont")}); copy it to '
                         f'reports/campaign-20261012/explore-ms/summaries-v2/{CHR22_RUN.format(h=h, s="simont")}/summary.json and rerun.')
    out.extend(['## Still missing', ''] + still + [''])

    out.extend(['## Gaps and questions', '', '### Cells absent or not checked (generated)', ''])
    gaps = missing_cells(ms, chr22, loc)
    if gaps:
        out.extend(f'- {m}' for m in gaps)
    else:
        out.append('- None among the plotted series.')
    out.extend(['', '### Naming', '',
                '- The ONT set in the Human472 MS and offset runs is human_ont_modern, HG002 Q28; panels use '
                '"Nanopore (HG002 Q28)".',
                '- Guccione ERR14066769 is an Illumina paired-end run with both mates entering the workload as '
                'interleaved reads (EXECUTION.md, October 16: 219,488 read pairs as 439,912 records); panels use '
                '"Skin metagenome (Illumina, ERR14066769)".',
                '- The short-read MS set is the ACGT-only 60 Mbase set (478,638 reads); the short-read offset set '
                'is a different, smaller ACGT-only set (49,152 reads, 6.14 Mbases).', '',
                '### Other observations', '',
                '- ms2-phi is out of the MS-only set and absent from every MS figure. Its cells remain in the '
                'summaries and in no plot.',
                '- Oracle placement: for ms2-dual the replay runs measured lower peaks than ordinary (about '
                '127 GiB against 184 GiB on Human472). The figures place replay ticks at the ordinary peak, as '
                'requested; `plotted-points.csv` keeps both values.',
                '- ms2-psi and ms2-dual locate from the same index file (human472.tms.full.seqids.orbitidx, '
                '362 GB), so their peaks are equal (about 293 GiB, 296 GiB with ftab). This also means the ms2-psi '
                'locate peak is far above its MS-only peak (97 GiB); a Psi-appropriate locate index is to be '
                'checked before building anything (October 25).',
                '- ropebwt3 -s6 (about 345 GiB) uses more memory than ms2 locate (about 293 GiB) and is still '
                '5 to 10 times slower in locations per second.',
                '- ftab in locate: on short reads ftab helps at every cap, but on PacBio and ONT at N = 100 the ftab '
                'variants are slower than without (for example PacBio Psi 3.40M against 3.75M locations/s). Likely '
                'node noise; only short reads are plotted with ftab.',
                '- Each Human472 read set ran on a different node (MS: bigmem06, 09, 10, 11; offsets: bigmem08, '
                '14, 16, 17). Earlier notes record node-to-node speed variation of a few percent, so cross-panel '
                'comparisons carry that noise.',
                '- MS v2 runs record three runs per cell (warmup plus two timed); the value plotted is the median '
                'of the two timed runs.',
                '- Fig Ops time is the query phase as the harness measures it (read, parse, query, format, '
                'write). Orbit\'s own query-only time is much lower on mtb152 (30.6 against 55.6 ns per base '
                'for ms2-lcp), so on mtb152 about 25 ns per base of every bar is I/O and parsing that no layout '
                'can remove; the relative gaps between layouts are larger in query-only terms. On HPRC2 full the '
                'same overhead is about 26 ns per base of 140.',
                '- Fig Ops peak memory on HPRC2 full: ms2-lcp 90 GiB, ms2-psi 102 GiB, ms2-phi 122 GiB, ms2-dual '
                '189 GiB (ordinary); replay runs load less (ms2-dual 132 GiB, ms2-phi 65 GiB). The figure does '
                'not draw memory; the Load times table and plotted-points.csv keep it.',
                '- Fig Ops ms2-phi uses a plum (#a3478f) that passes the palette checker with the other three '
                'layout colors but is too close to the Others purple and ropebwt3 red to share a figure with '
                'them; no current figure does.',
                '- Figure colors: layout colors, ropebwt3 red and the Others purple were validated as one set with '
                'the palette checker (no normal-vision failure). Two colors (psi aqua, dual dark yellow) sit just '
                'under 3:1 contrast on white, and the closest color-vision-deficiency pair (dual against ropebwt3 '
                'red) is in the warning band, so identity also rides on markers, direct labels and the legend. '
                'Grayscale printing will not separate the layouts.',
                ''])
    (out_dir / 'CAPTIONS.md').write_text('\n'.join(out))


def write_csv(rows, out_dir):
    if not rows:
        return
    with open(out_dir / 'plotted-points.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--ms-dir', default=DEFAULT_MS)
    ap.add_argument('--chr22-dir', default=DEFAULT_CHR22)
    ap.add_argument('--locate-dir', default=DEFAULT_LOCATE)
    ap.add_argument('--ops-dir', default=DEFAULT_OPS)
    ap.add_argument('--ops2-reads', choices=sorted(OPS2_RUNS), default=OPS2_DEFAULT_READS,
                    help='simulated read type drawn in fig-ops2')
    ap.add_argument('--out', default=DEFAULT_OUT)
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    ms = load_ms(args.ms_dir)
    chr22 = load_chr22(args.chr22_dir)
    loc = load_locate(args.locate_dir)
    ops = load_ops(args.ops_dir)
    ops2 = load_ops(args.ops_dir, OPS2_RUNS[args.ops2_reads])
    rec = Recorder()
    ms_grid_figure(ms, 'fig-ms', out_dir, rec, replay=False)
    ms_grid_figure(ms, 'fig-ms-oracle', out_dir, rec, replay=True)
    if 'human_short' in ms:
        ms_ftab_figure(ms, out_dir, rec)
    chr22_drawn = [h for h in CHR22_REFS if chr22_figure(chr22, h, out_dir, rec)]
    for cap in CAPS:
        locate_grid_figure(loc, cap, out_dir, rec)
    if 'human_short' in loc:
        locate_ftab_figure(loc, out_dir, rec)
    for cap in LORACLE_CAPS:
        locate_oracle_figure(loc, cap, out_dir, rec)
    ops_figure(ops, out_dir, rec)
    ops_figure(ops2, out_dir, rec, stem='fig-ops2')
    write_csv(rec.rows, out_dir)
    write_captions(ms, chr22, loc, ops, rec.rows, chr22_drawn, out_dir, ops2, args.ops2_reads)
    for h in CHR22_REFS:
        gone = [s for s in CHR22_SETS if (h, s) not in chr22]
        if gone:
            print(f'fig-ms22-{h}: missing {", ".join(gone)}')
    print(f'wrote {out_dir}')


if __name__ == '__main__':
    main()
