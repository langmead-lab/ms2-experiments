"""Check every experimental number the paper quotes against the committed results.

Usage:
  python -m ms2exp.figures.paper_numbers [--figures figures] [--results results]
                                         [--read-sets data/read-sets.json]
                                         [--values figures/paper-values.json] [--rep-tsv FILE]
                                         [--md OUT.md] [--table OUT.tex] [--paper DIR] [--strict]

figures/paper-values.json lists each number as the paper prints it, by id. This script computes
every id from the committed records, rounds it half up to the precision the paper prints, and
compares:

  figures/draft/plotted-points.csv      plotted speeds, peaks, loads and queries
  figures/summaries/GROUP/summary.json  workloads, count cells, per-run loads and hosts
  results/MANIFEST/.../summary.json     the transparent huge pages (THP) comparison and thread scaling
  data/read-sets.json                   read counts and length distributions

An entry with expected_change is a known disagreement between the paper and the data: the data must
then agree with expected_change.to (or still disagree with the paper, when to is null), and the
entry's edit is the wording the paper needs. Without one, the data must agree with the paper.
--strict exits nonzero when any check fails, so the text cannot drift from the data unnoticed.

--paper DIR checks the other direction: every entry's number (or its text field, for an entry whose
claim is worded rather than numeric) must still appear in the paper's LaTeX source in DIR, outside
comments, so an entry cannot keep agreeing with text the paper no longer has. With --strict, a
missing entry fails the run. The paper commit checked is printed, with a note when it is not the one
the values file names.

--md writes the report (figures/draft/paper-numbers.md from the figures driver), --table Table A2
(index load times) as a LaTeX tabular. Fig rep's means come from runs outside the freeze
(results/README.md) and are checked only with --rep-tsv, a fig-rep-batches.tsv of those runs.
"""
import argparse
import csv
from collections import defaultdict
from decimal import Decimal, ROUND_HALF_UP
import glob
import json
from pathlib import Path
import re
import statistics
import subprocess
import sys

REPO = Path(__file__).resolve().parents[2]
FREEZE = 'freeze-2885ee8-00e665de'
THP_OFF = 'freeze-2885ee8-3c8c8967'
THREADS = 'threads-2885ee8-da1041e8'
HUMAN = ['human_short', 'human_pacbio', 'human_ont_modern']
COHORTS = HUMAN + ['guccione_skin']
SIMS = ['simshort', 'simpacbio', 'simont']
LEGACY = ['MONI', 'PHONI', 'AugPHONI']
MS_GROUP = '{}-ms-acgt60mb'
LOC_GROUPS = ['human_short-locate-6mb', 'human_pacbio-locate-20mb', 'human_ont_modern-locate-10mb',
              'guccione_skin-locate-500kb']
OPS_GROUPS = {'mtb152': 'fig_ops-mtb152-heldout16-p100k-s600-v1',
              'human472': 'fig_ops-human472-apr008-hap1-p100k-s600-v1',
              'mtb152 reads': 'fig_ops-mtb152-pacbio-60mb-v2',
              'human472 reads': 'fig_ops-human472-apr008-pacbio-60mb-v2'}
FLOW_CELL_MBASES = 2675 * 420 / 1e6  # channels per PromethION flow cell times bases per second per channel


class Value:
    """A computed number or set of numbers and how the paper states it.

kind: value (one number), range (the minimum and maximum of numbers), atleast (the minimum, which
the paper's bound must not exceed by a whole unit), within (the maximum, which must round to at most
the paper's bound), about (one number within 5% of the paper's).
"""

    def __init__(self, kind, numbers):
        self.kind = kind
        self.numbers = numbers if isinstance(numbers, (list, tuple)) else [numbers]


def decimals(text):
    text = text.replace(',', '')
    return len(text.split('.')[1]) if '.' in text else 0


def rnd(x, nd):
    """x rounded half up to nd decimals, from its shortest decimal form."""
    return Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP)


def number(text):
    return float(text.replace(',', ''))


def show(v, text):
    """The computed value written at the precision of text."""
    if v.kind == 'range':
        ends = text.split('--') if '--' in text else [text, text]
        return f'{rnd(min(v.numbers), decimals(ends[0])):,} to {rnd(max(v.numbers), decimals(ends[1])):,}'
    if v.kind == 'atleast':
        return f'min {min(v.numbers):.2f}'
    if v.kind == 'within':
        return f'max {max(v.numbers):.2f}'
    if v.kind == 'about':
        return f'{v.numbers[0]:.2f}'
    return f'{rnd(v.numbers[0], decimals(text)):,}'


def agrees(v, text):
    """Whether the value agrees with text at text's precision; None when text is not a number."""
    try:
        if v.kind == 'range':
            lo_t, hi_t = text.split('--') if '--' in text else (text, text)
            return (rnd(min(v.numbers), decimals(lo_t)) == rnd(number(lo_t), decimals(lo_t)) and
                    rnd(max(v.numbers), decimals(hi_t)) == rnd(number(hi_t), decimals(hi_t)))
        bound = number(text)
    except ValueError:
        return None
    if v.kind == 'atleast':
        return bound <= min(v.numbers) < bound + 1
    if v.kind == 'within':
        return rnd(max(v.numbers), decimals(text)) <= rnd(bound, decimals(text))
    if v.kind == 'about':
        return abs(v.numbers[0] - bound) / bound < 0.05
    return rnd(v.numbers[0], decimals(text)) == rnd(bound, decimals(text))


class Data:
    def __init__(self, figures, results, read_sets):
        self.results = Path(results)
        self.pp = {}
        with open(Path(figures)/'draft/plotted-points.csv') as fh:
            for r in csv.DictReader(fh):
                self.pp[(r['figure'], r['read_set'], r['series'])] = r
        self.groups = {p.parent.name: json.loads(p.read_text()) for p in (Path(figures)/'summaries').glob('*/summary.json')}
        self.read_sets = json.loads(Path(read_sets).read_text())['sets']

    def speed(self, fig, rs, series):
        return float(self.pp[(fig, rs, series)]['speed'])

    def mem(self, fig, rs, series):
        return float(self.pp[(fig, rs, series)]['peak_gib_measured'])

    def cell(self, group, cid):
        return next(c for c in self.groups[group]['cells'] if c['id'] == cid)

    def workload(self, group):
        return self.groups[group]['config']['workload']

    def read_set(self, group):
        w = self.workload(group)
        return self.read_sets[f"{w['reference']}/{w['read_set']}"]

    def run_cells(self, manifest):
        out = {}
        for p in glob.glob(str(self.results/manifest/'*/cells/*/summary.json')):
            d = json.loads(Path(p).read_text())
            out[(d['job'], d['cell'])] = d
        return out


def ops_values(D, v):
    rep = {k: D.cell(g + '-counts', 'ms-ordinary')['ops_per_base']['line']['repositions/base'] for k, g in OPS_GROUPS.items()}
    for k, g in OPS_GROUPS.items():
        # Recoveries belong to the query and index, so every count cell of a group must agree.
        if len({round(c['ops_per_base']['line']['repositions/base'], 9) for c in D.groups[g + '-counts']['cells']}) != 1:
            raise ValueError(f'{g}-counts: the count cells disagree on repositions per base')
    v['ops.mtb152_genome_mbases'] = Value('value', D.workload(OPS_GROUPS['mtb152'])['bases'] / 1e6)
    v['ops.hprc2_genome_mbases'] = Value('value', D.workload(OPS_GROUPS['human472'])['bases'] / 1e6)
    v['ops.pieces'] = Value('range', [D.workload(OPS_GROUPS[k])['reads'] for k in ('mtb152', 'human472')])
    v['ops.hifi_mbases'] = Value('range', [D.workload(OPS_GROUPS[k])['bases'] / 1e6 for k in ('mtb152 reads', 'human472 reads')])
    v['ops.hifi_median_kbases'] = Value('range', [D.read_set(OPS_GROUPS[k])['length']['median'] / 1e3
                                                  for k in ('mtb152 reads', 'human472 reads')])
    v['ops.rec_mtb152_genomes'] = Value('value', rep['mtb152'])
    v['ops.rec_hprc2_genomes'] = Value('value', rep['human472'])
    v['ops.rec_mtb152_hifi'] = Value('value', rep['mtb152 reads'])
    v['ops.rec_hprc2_hifi'] = Value('value', rep['human472 reads'])
    v['ops.rec_hifi_over_genomes'] = Value('range', [rep['mtb152 reads'] / rep['mtb152'],
                                                     rep['human472 reads'] / rep['human472']])
    ns = lambda fig, ref, s: D.speed(fig, ref, s)
    gap = lambda fig, ref, s: 100 * (ns(fig, ref, s) / ns(fig, ref, s + ' oracle') - 1)
    v['ops.genomes_min_dual_gap_hprc2'] = Value('within', [gap('fig-ops', 'human472', s) for s in ('ms2-min', 'ms2-dual')])
    for s in ('min', 'dual', 'psi', 'phi'):
        v[f'ops.hprc2_{s}_ns'] = Value('value', ns('fig-ops', 'human472', f'ms2-{s}'))
        v[f'ops.mtb152_{s}_ns'] = Value('value', ns('fig-ops', 'mtb152', f'ms2-{s}'))
    for s in ('min', 'dual'):
        v[f'ops.hprc2_{s}_oracle_ns'] = Value('value', ns('fig-ops', 'human472', f'ms2-{s} oracle'))
    for s in ('dual', 'psi', 'phi'):
        v[f'ops.hprc2_{s}_ops'] = Value('value', ns('fig-ops', 'human472', f'ms2-{s} ops'))
    v['ops.hifi_min_gap'] = Value('within', [gap('fig-ops2', r, 'ms2-min') for r in ('mtb152', 'human472')])
    v['ops.hifi_dual_gap'] = Value('range', [gap('fig-ops2', r, 'ms2-dual') for r in ('mtb152', 'human472')])
    v['ops.hifi_hprc2_phi_over_oracle'] = Value('value', 1 + gap('fig-ops2', 'human472', 'ms2-phi') / 100)


def chr22_values(D, v):
    f100, f10 = 'fig-ms22-100', 'fig-ms22-10'
    mn = [D.speed(f100, s, 'ms2-min') for s in SIMS]
    best = [max(D.speed(f100, s, t) for t in LEGACY) for s in SIMS]
    phi = [D.speed(f100, s, 'ms2-phi') for s in SIMS]
    v['chr22.min_over_phoni'] = Value('range', [a / D.speed(f100, s, 'PHONI') for a, s in zip(mn, SIMS)])
    v['chr22.min_mbases'] = Value('range', mn)
    v['chr22.min_over_best_legacy'] = Value('range', [a / b for a, b in zip(mn, best)])
    v['chr22.min_gib'] = Value('range', [D.mem(f100, s, 'ms2-min') for s in SIMS])
    v['chr22.legacy_gib'] = Value('range', [D.mem(f100, s, t) for s in SIMS for t in LEGACY])
    v['chr22.phi_mbases'] = Value('range', phi)
    v['chr22.phi_over_best_legacy'] = Value('range', [a / b for a, b in zip(phi, best)])
    v['chr22.rope_gib'] = Value('value', max(D.mem(f100, s, 'ropebwt3 -l1') for s in SIMS))
    v['chr22.min_over_rope'] = Value('range', [a / D.speed(f100, s, 'ropebwt3 -l1') for a, s in zip(mn, SIMS)])
    v['chr22.h10_change_min_psi_dual'] = Value('within', [100 * abs(D.speed(f10, s, l) / D.speed(f100, s, l) - 1)
                                                           for s in SIMS for l in ('ms2-min', 'ms2-psi', 'ms2-dual')])
    v['chr22.phi_h10_over_h100'] = Value('range', [D.speed(f10, s, 'ms2-phi') / D.speed(f100, s, 'ms2-phi') for s in SIMS])
    v['chr22.phi_h10_mbases'] = Value('range', [D.speed(f10, s, 'ms2-phi') for s in SIMS])
    v['chr22.reads_mbases'] = Value('range', [D.workload(f'chr22-h100-ms-{s}')['bases'] / 1e6 for s in SIMS])
    mn10 = [D.speed(f10, s, 'ms2-min') for s in SIMS]
    v['app10.min_mbases'] = Value('range', mn10)
    v['app10.min_change'] = Value('within', [100 * abs(a / b - 1) for a, b in zip(mn10, mn)])
    v['app10.legacy_mbases'] = Value('range', [D.speed(f10, s, t) for s in SIMS for t in LEGACY])
    v['app10.min_gib'] = Value('range', [D.mem(f10, s, 'ms2-min') for s in SIMS])
    v['app10.legacy_gib'] = Value('range', [D.mem(f10, s, t) for s in SIMS for t in LEGACY])


def hprc2_values(D, v):
    f = 'fig-ms-oracle'
    s = {(c, m): D.speed(f, c, m) for (fig, c, m) in D.pp if fig == f}
    rope = {c: s[(c, 'ropebwt3 -l1')] for c in COHORTS}
    for mode in ('min', 'dual', 'psi'):
        v[f'ms.{mode}_over_rope_human'] = Value('atleast', [s[(c, f'ms2-{mode}')] / rope[c] for c in HUMAN])
    v['ms.rope_size_pct'] = Value('about', 100 * statistics.median(D.mem(f, c, 'ropebwt3 -l1') for c in HUMAN) /
                                  statistics.median(D.mem(f, c, 'ms2-min') for c in HUMAN))
    v['ms.oracle_mbases'] = Value('range', [x for (c, m), x in s.items() if m.endswith(' oracle')])
    v['ms.min_oracle_slowdown_human'] = Value('range', [s[(c, 'ms2-min oracle')] / s[(c, 'ms2-min')] for c in HUMAN])
    v['ms.psi_dual_oracle_slowdown_human'] = Value('range', [s[(c, m + ' oracle')] / s[(c, m)] for c in HUMAN
                                                             for m in ('ms2-psi', 'ms2-dual')])
    sk = 'guccione_skin'
    v['ms.skin_min'] = Value('value', s[(sk, 'ms2-min')])
    v['ms.skin_min_oracle'] = Value('value', s[(sk, 'ms2-min oracle')])
    v['ms.skin_psi'] = Value('value', s[(sk, 'ms2-psi')])
    v['ms.skin_dual'] = Value('value', s[(sk, 'ms2-dual')])
    v['ms.skin_psi_dual_oracle'] = Value('about', statistics.mean([s[(sk, 'ms2-psi oracle')], s[(sk, 'ms2-dual oracle')]]))
    v['ms.phi_mbases'] = Value('range', [s[(c, 'ms2-phi')] for c in COHORTS])
    v['ms.phi_pct_of_psi'] = Value('range', [100 * s[(c, 'ms2-phi')] / s[(c, 'ms2-psi')] for c in COHORTS])
    v['ms.min_over_rope_all'] = Value('range', [s[(c, 'ms2-min')] / rope[c] for c in COHORTS])
    v['ms.rope_l1_mbases'] = Value('range', list(rope.values()))
    v['ms.phi_over_rope'] = Value('range', [s[(c, 'ms2-phi')] / rope[c] for c in COHORTS])
    v['ms.min_ont_mbases'] = Value('value', s[('human_ont_modern', 'ms2-min')])
    v['ms.flow_cells'] = Value('value', s[('human_ont_modern', 'ms2-min')] / FLOW_CELL_MBASES)
    for c in COHORTS:
        w = D.workload(MS_GROUP.format(c))
        v[f'reads.{c}_reads'] = Value('value', w['reads'])
        v[f'reads.{c}_mbases'] = Value('value', w['bases'] / 1e6)
    v['reads.all_mbases'] = Value('range', [D.workload(MS_GROUP.format(c))['bases'] / 1e6 for c in COHORTS])
    pb = D.read_set(MS_GROUP.format('human_pacbio'))['length']
    ont = D.read_set(MS_GROUP.format('human_ont_modern'))['length']
    v['reads.pacbio_median_kbases'] = Value('value', pb['median'] / 1e3)
    v['reads.ont_min_bp'] = Value('value', ont['min'])
    v['reads.ont_max_kbases'] = Value('value', ont['max'] / 1e3)
    v['reads.ont_median_kbases'] = Value('value', ont['median'] / 1e3)
    v['reads.ont_n50_kbases'] = Value('value', ont['n50'] / 1e3)


def locate_values(D, v):
    loc = lambda cap, c, m: D.speed(f'fig-l{cap}', c, m)
    ms2 = ('ms2-psi', 'ms2-dual')
    for g, key in zip(LOC_GROUPS, COHORTS):
        v[f'loc.{key}_reads'] = Value('value', D.workload(g)['reads'])
        v[f'loc.{key}_mbases'] = Value('value', D.workload(g)['bases'] / 1e6)
    v['loc.ms2_gib'] = Value('value', statistics.median(D.mem('fig-l100', c, m) for c in COHORTS
                                                        for m in ('ms2-psi', 'ms2-phi', 'ms2-dual')))
    v['loc.n100_psi_dual_mloc'] = Value('range', [loc(100, c, m) / 1e6 for c in COHORTS for m in ms2])
    v['loc.n100_over_s8'] = Value('range', [loc(100, c, m) / loc(100, c, 'ropebwt3 -s8') for c in COHORTS for m in ms2])
    v['loc.n100_over_s6'] = Value('range', [loc(100, c, m) / loc(100, c, 'ropebwt3 -s6') for c in COHORTS for m in ms2])
    for shift in (8, 6, 11):
        v[f'loc.s{shift}_gib'] = Value('value', statistics.median(D.mem('fig-l100', c, f'ropebwt3 -s{shift}') for c in COHORTS))
    v['loc.n100_phi_mloc'] = Value('range', [loc(100, c, 'ms2-phi') / 1e6 for c in COHORTS])
    v['loc.n100_short_phi_over_dual'] = Value('value', loc(100, 'human_short', 'ms2-phi') / loc(100, 'human_short', 'ms2-dual'))
    v['loc.n100_other_phi_over_dual'] = Value('range', [loc(100, c, 'ms2-phi') / loc(100, c, 'ms2-dual') for c in COHORTS[1:]])
    v['loc.n100_s11_kloc'] = Value('range', [loc(100, c, 'ropebwt3 -s11') / 1e3 for c in COHORTS])
    v['appn.over_s6'] = Value('range', [loc(cap, c, m) / loc(cap, c, 'ropebwt3 -s6') for cap in (1, 10)
                                        for c in COHORTS for m in ms2])
    v['appn.n1_kloc'] = Value('range', [loc(1, c, m) / 1e3 for c in COHORTS for m in ms2])
    v['appn.n10_mloc'] = Value('range', [loc(10, c, m) / 1e6 for c in COHORTS for m in ms2])
    v['appn.n1_short_phi_over_dual'] = Value('value', loc(1, 'human_short', 'ms2-phi') / loc(1, 'human_short', 'ms2-dual'))
    v['appn.n10_short_phi_over_dual'] = Value('value', loc(10, 'human_short', 'ms2-phi') / loc(10, 'human_short', 'ms2-dual'))
    v['appn.other_phi_over_dual'] = Value('range', [loc(cap, c, 'ms2-phi') / loc(cap, c, 'ms2-dual') for cap in (1, 10)
                                                    for c in COHORTS[1:]])


# Table A2: (row name, use, group kind, cell ids). Load is the median over every run of the cells
# behind Figs ms-oracle (ordinary, no ftab, ms2-dual on the full index with sequence names) and
# locate100 (N = 100); the range is over per-node medians; resident is the median peak of the cells.
LOAD_ROWS = [('ms2-min', 'matching statistics', 'ms', ['ms-ordinary']),
             ('ms2-psi', 'matching statistics', 'ms', ['psi-ordinary']),
             ('ms2-phi', 'matching statistics', 'ms', ['phi-ordinary']),
             ('ms2-dual', 'matching statistics', 'ms', ['dual-ordinary']),
             ('ms2-psi, ms2-phi, ms2-dual', 'locating', 'loc', ['cap100-psi-ordinary', 'cap100-phi-ordinary',
                                                                'cap100-dual-ordinary']),
             ('ropebwt3', 'matching statistics', 'ms', ['rope-l1']),
             (r'ropebwt3 \texttt{-s11}', 'locating', 'loc', ['cap100-rope-s11']),
             (r'ropebwt3 \texttt{-s8}', 'locating', 'loc', ['cap100-rope-s8']),
             (r'ropebwt3 \texttt{-s6}', 'locating', 'loc', ['cap100-rope-s6'])]
LOAD_IDS = ['min', 'psi', 'phi', 'dual', 'ms2_locate', 'rope', 'rope_s11', 'rope_s8', 'rope_s6']


def load_table(D):
    groups = {'ms': [MS_GROUP.format(c) for c in COHORTS], 'loc': LOC_GROUPS}
    out = []
    for name, use, kind, ids in LOAD_ROWS:
        loads, peaks = [], []
        for g in groups[kind]:
            for c in D.groups[g]['cells']:
                if c['id'] in ids:
                    loads += c['load_seconds']
                    peaks.append(c['peak_rss_kib'] / 2 ** 20)
        median, resident = statistics.median(loads), statistics.median(peaks)
        out.append(dict(name=name, use=use, resident=resident, load=median,
                        rate=resident / median))
    return out


def load_table_tex(table):
    lines = ['% Written by python -m ms2exp.figures.paper_numbers from figures/summaries. Load: median over every',
             '% run of the cells behind Figs ms-oracle (ordinary, no ftab; ms2-dual on the full index with',
             '% sequence names; ropebwt3 -l1) and locate100 (N=100).',
             '% Resident: median peak RSS of those cells. GiB/s: resident over median load.',
             r'\begin{tabular}{llrrr}', r'\hline', r'Index & Used for & Resident (GiB) & Load (s) & GiB/s \\', r'\hline']
    for r in table:
        lines.append(f"{r['name']} & {r['use']} & {rnd(r['resident'], 0)} & {rnd(r['load'], 1)} & {rnd(r['rate'], 1)} \\\\")
    return '\n'.join(lines + [r'\hline', r'\end{tabular}', ''])


def load_values(D, v, table):
    for key, r in zip(LOAD_IDS, table):
        v[f'load.{key}_resident'] = Value('value', r['resident'])
        v[f'load.{key}_seconds'] = Value('value', r['load'])
        v[f'load.{key}_rate'] = Value('value', r['rate'])
    v['load.ms2_rates'] = Value('range', [r['rate'] for r in table if r['name'].startswith('ms2')])
    v['load.rope_rates'] = Value('about', statistics.median(r['rate'] for r in table if r['name'].startswith('rope')))
    ratio = {c['id']: c['median_load_seconds'] / c['median_query_seconds']
             for c in D.groups[MS_GROUP.format('human_short')]['cells']}
    v['load.short_load_over_query'] = Value('range', [ratio['ms-ordinary'], ratio['psi-ordinary'], ratio['dual-ordinary']])
    v['load.short_load_over_query_phi'] = Value('value', ratio['phi-ordinary'])


def thp_comparison(D):
    """Cells in both freeze runs with three timed runs in each, and the cells the rerun copied."""
    off, on = D.run_cells(THP_OFF), D.run_cells(FREEZE)
    common = [k for k in off if k in on]
    timed = [k for k in common if len(off[k].get('query_seconds') or []) == 3 and len(on[k].get('query_seconds') or []) == 3]
    copied = [k for k in common if k not in timed]
    if any(off[k]['query_seconds'] != on[k]['query_seconds'] for k in copied):
        raise ValueError('a cell with fewer than three runs differs between the runs, so it was not copied')
    return off, on, common, timed, copied


def thp_values(D, v, notes):
    off, on, common, timed, copied = thp_comparison(D)
    by, mem = defaultdict(list), []
    for k in timed:
        a, b = off[k], on[k]
        if a['reference'] == 'human472':
            by[(a['tool'], 'locate' if 'locate' in a['group'] else 'ms')].append(a['median_query_seconds'] / b['median_query_seconds'])
            mem.append(100 * (b['peak_rss_kib'] / a['peak_rss_kib'] - 1))
    v['thp.ms2_range'] = Value('range', by[('ms2', 'ms')] + by[('ms2', 'locate')])
    v['thp.ms2_ms_median'] = Value('value', statistics.median(by[('ms2', 'ms')]))
    v['thp.ms2_locate_median'] = Value('value', statistics.median(by[('ms2', 'locate')]))
    v['thp.rope_median'] = Value('about', statistics.median(by[('ropebwt3', 'ms')] + by[('ropebwt3', 'locate')]))
    v['thp.peak_change_pct'] = Value('within', [abs(x) for x in mem])
    rope = by[('ropebwt3', 'ms')] + by[('ropebwt3', 'locate')]
    notes.append(f'THP comparison: {len(common)} cells in both runs, {len(timed)} timed three times in each '
                 f'({len(copied)} single-run count and Fig rep cells were copied into the rerun, identical); '
                 f'{sum(map(len, by.values()))} of the timed cells are on HPRC2 (ms2 MS {len(by[("ms2", "ms")])}, '
                 f'ms2 locate {len(by[("ms2", "locate")])}, ropebwt3 MS {len(by[("ropebwt3", "ms")])}, ropebwt3 '
                 f'locate {len(by[("ropebwt3", "locate")])}); ropebwt3 speedups {min(rope):.2f} to {max(rope):.2f}.')


def thread_values(D, v):
    cells = [d for d in D.run_cells(THREADS).values() if d['group'].startswith('threads-human_')]
    sp = {(d['read_set'], d.get('layout') or 'rope', d['threads']): d['mbases_per_second'] for d in cells}
    sets = sorted({k[0] for k in sp})
    w = {s: next(d for d in cells if d['read_set'] == s) for s in sets}
    for s in sets:
        key = s.split('-ms-')[0]
        v[f'threads.{key}_reads'] = Value('value', w[s]['reads'])
        v[f'threads.{key}_mbases'] = Value('value', w[s]['bases'] / 1e6)
    modes = ['ms', 'psi', 'phi', 'dual']
    v['threads.ms2_32_over_1'] = Value('range', [sp[(s, m, 32)] / sp[(s, m, 1)] for s in sets for m in modes])
    v['threads.rope_32_over_1'] = Value('range', [sp[(s, 'rope', 32)] / sp[(s, 'rope', 1)] for s in sets])
    v['threads.min_32_mbases'] = Value('range', [sp[(s, 'ms', 32)] for s in sets])
    v['threads.phi_rope_64_over_32'] = Value('range', [sp[(s, m, 64)] / sp[(s, m, 32)] for s in sets for m in ('phi', 'rope')])
    v['threads.min_dual_64_over_32'] = Value('within', [sp[(s, m, 64)] / sp[(s, m, 32)] for s in sets for m in ('ms', 'dual')])


def rep_values(tsv, v):
    means = defaultdict(dict)
    with open(tsv) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['mode'] in ('psi', 'phi', 'dual'):
                means[(r['series'], r['mode'])][r['run']] = float(r['mean'])
    a, b = (lambda m: means[('a2', m)]), (lambda m: means[('b2', m)])
    v['rep.a_psi_exact'] = Value('value', a('psi')['exact'])
    v['rep.a_psi_014'] = Value('value', a('psi')['e0.0014'])
    v['rep.a_psi_floor'] = Value('value', a('psi')['e0.020'])
    v['rep.a_phi'] = Value('range', list(a('phi').values()))
    v['rep.b_phi'] = Value('range', [b('phi')['h5'], b('phi')['h464']])
    v['rep.b_psi'] = Value('range', list(b('psi').values()))
    v['rep.b_dual'] = Value('range', list(b('dual').values()))


def compute(figures, results, read_sets, rep_tsv=None):
    D = Data(figures, results, read_sets)
    v, notes = {}, []
    ops_values(D, v)
    chr22_values(D, v)
    hprc2_values(D, v)
    locate_values(D, v)
    table = load_table(D)
    load_values(D, v, table)
    thp_values(D, v, notes)
    thread_values(D, v)
    if rep_tsv:
        rep_values(rep_tsv, v)
    return v, table, notes


def check(values, spec):
    """One row per entry: (entry, computed text, status). status is agrees, edit, unchecked or FAIL."""
    rows = []
    for e in spec['numbers']:
        val = values.get(e['id'])
        if val is None:
            if not e.get('optional'):
                raise KeyError(f"no computation for {e['id']}")
            rows.append((e, 'data outside the repository', 'unchecked'))
            continue
        change = e.get('expected_change')
        if change is None:
            rows.append((e, show(val, e['paper']), 'agrees' if agrees(val, e['paper']) else 'FAIL'))
            continue
        target = change.get('to')
        still_wrong = agrees(val, e['paper']) is not True
        ok = still_wrong and (target is None or agrees(val, target))
        rows.append((e, show(val, target or e['paper']), 'edit' if ok else 'FAIL'))
    unknown = set(values) - {e['id'] for e in spec['numbers']}
    if unknown:
        raise KeyError(f'computed ids missing from the values file: {sorted(unknown)}')
    return rows


def paper_text(paper_dir):
    """The paper's LaTeX outside comments, with thin spaces, ties, {,} and \\% made plain."""
    paper_dir = Path(paper_dir)
    parts = []
    for f in sorted(paper_dir.glob('*.tex')):
        parts.append(re.sub(r'(?m)(?<!\\)%.*$', '', f.read_text()))
    text = '\n'.join(parts)
    for a, b in (('{,}', ','), ('\\,', ' '), ('~', ' '), ('\\%', '%')):
        text = text.replace(a, b)
    return re.sub(r'\s+', ' ', text)


def text_forms(entry):
    if entry.get('text'):
        return [entry['text']]
    v = entry['paper']
    if '--' in v:
        a, b = v.split('--')
        return [f'{a} to {b}', f'{a}--{b}', f'{a} and {b}', f'{a}-{b}']
    return [v]


def missing_from_text(spec, text):
    """The entries whose number or text field does not appear in the paper text."""
    def found(form):
        return re.search(r'(?<![\d.])' + re.escape(form) + r'(?![\d])', text) is not None
    return [e for e in spec['numbers'] if not any(found(f) for f in text_forms(e))]


def paper_commit(paper_dir):
    try:
        out = subprocess.run(['git', '-C', str(paper_dir), 'rev-parse', '--short', 'HEAD'],
                             capture_output=True, text=True, check=True)
        dirty = subprocess.run(['git', '-C', str(paper_dir), 'status', '--porcelain', '--untracked-files=no'],
                               capture_output=True, text=True, check=True).stdout.strip()
        return out.stdout.strip() + (' with uncommitted edits' if dirty else '')
    except (OSError, subprocess.CalledProcessError):
        return None


def paper(text):
    return text.replace('--', ' to ')


def report(rows, spec, notes):
    out = ['# Paper numbers', '',
           'Written by `python -m ms2exp.figures.paper_numbers` from the committed results. Each number the paper '
           f"quotes ({spec['paper']}) is computed from the run records and rounded to the precision the paper prints. "
           'Status: agrees; edit (a known disagreement, with the paper edit below); FAIL (the paper and the data '
           'disagree, or a recorded edit no longer applies); unchecked (the data are outside this repository).', '']
    section = None
    for e, got, status in rows:
        if e['where'] != section:
            section = e['where']
            out += ['', f'## {section}', '', '| Number | Paper | Computed | Status |', '|---|---|---|---|']
        shown = paper(e['paper']) + (f" (to be {paper(e['expected_change']['to'])})"
                                     if e.get('expected_change', {}).get('to') else '')
        out.append(f"| {e['label']} | {shown} | {got} | {status} |")
    out += ['', '## Not checked here', '']
    out += [f"- {u['where']}: {u['label']}, \"{u['paper']}\". {u['reason']}" for u in spec['unchecked']]
    out += ['', '## Notes', ''] + [f'- {n}' for n in notes]
    edits = [e for e, _, status in rows if status in ('edit', 'unchecked') and e.get('expected_change')]
    out += ['', '## Paper edits outstanding', '']
    items = [(e['where'], e['expected_change']['edit']) for e in edits] + [(t['where'], t['edit'])
                                                                           for t in spec.get('text_edits', [])]
    out += [f'{i}. {where}: {text}' for i, (where, text) in enumerate(items, 1)] or ['None.']
    counts = defaultdict(int)
    for _, _, status in rows:
        counts[status] += 1
    out += ['', 'Totals: ' + ', '.join(f'{counts[s]} {s}' for s in ('agrees', 'edit', 'unchecked', 'FAIL')) + '.', '']
    return '\n'.join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--figures', default=REPO/'figures')
    ap.add_argument('--results', default=REPO/'results')
    ap.add_argument('--read-sets', default=REPO/'data/read-sets.json')
    ap.add_argument('--values', default=REPO/'figures/paper-values.json')
    ap.add_argument('--rep-tsv')
    ap.add_argument('--md')
    ap.add_argument('--table')
    ap.add_argument('--paper', help="the paper's LaTeX source directory, to check that each number is still in it")
    ap.add_argument('--strict', action='store_true')
    args = ap.parse_args()
    spec = json.loads(Path(args.values).read_text())
    values, table, notes = compute(args.figures, args.results, args.read_sets, args.rep_tsv)
    rows = check(values, spec)
    text = report(rows, spec, notes)
    if args.md:
        Path(args.md).write_text(text)
    if args.table:
        Path(args.table).parent.mkdir(parents=True, exist_ok=True)
        Path(args.table).write_text(load_table_tex(table))
    failed = [(e, got) for e, got, status in rows if status == 'FAIL']
    for e, got in failed:
        print(f"FAIL {e['id']}: paper {e['paper']}, computed {got}", file=sys.stderr)
    print(text.rsplit('\n', 2)[-2])
    missing = []
    if args.paper:
        missing = missing_from_text(spec, paper_text(args.paper))
        commit = paper_commit(args.paper)
        pinned = spec['paper'].split()[-1]
        print(f"Paper text at {commit or 'an unknown commit'}: {len(spec['numbers']) - len(missing)} entries found, "
              f"{len(missing)} missing.")
        if commit and not commit.startswith(pinned):
            print(f'Note: the values file names {pinned}; update its paper field after checking.')
        for e in missing:
            print(f"MISSING {e['id']}: {' / '.join(text_forms(e))} not in the paper text", file=sys.stderr)
    return 1 if (failed or missing) and args.strict else 0


if __name__ == '__main__':
    sys.exit(main())
