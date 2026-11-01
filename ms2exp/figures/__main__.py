"""Regenerate the paper's figure data, draft figures and TikZ figures from the committed results.

Usage:
  python -m ms2exp.figures [--results results] [--out figures] [--manifest freeze-2885ee8-00e665de]
                           [--no-pdf]

Steps, each written under OUT:
  summaries/       ms2exp.collect --no-compare --platform aws on the manifest's results: one
                   GROUP/summary.json per cell group, and report.json
  draft/           paper_figures.py: matplotlib drafts (PNG, SVG, PDF), CAPTIONS.md and
                   plotted-points.csv, every plotted value with its load time and node
  tikz/            the TikZ sources (copied from figures/tikz when OUT is elsewhere)
  tikz/data/       paper_figures_tikz.py and gap_annotations.py: the data the TikZ figures read;
                   fig-rep/ from fig_rep_tikz.py --design 4 (the freeze's A1 runs, box statistics only)
  tex/             structure_bits.py's table of field widths (structure_bits.tex); the other .tex
                   files there are static snippets
  pdf/             the TikZ figures compiled with latexmk under the names the paper includes,
                   unless --no-pdf or latexmk is missing

paper_figures.py looks summaries up by run name, so a staging directory maps those names onto the
freeze groups: the MS, chr22 and locate groups by link, and each Fig ops run by merging a timed
group with its -counts group (ns per base from the timed cells, operations per base from the
counting cells).
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from ms2exp.figures.paper_figures import LOCATE_RUNS, OPS2_RUNS, OPS_RUNS

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
# The TikZ sources, copied from the paper's figs/tikz; they read their data from data/ beside them.
SOURCES = REPO/'figures/tikz'
TIKZ_SOURCES = ('figures-style.tex', 'figures-tweaks.tex', 'figures-gap.tex', 'fig-ops.tex', 'fig-ops-row.tex',
                'fig-ms22-100.tex', 'fig-ms22-10.tex', 'fig-ms-oracle.tex', 'fig-l1.tex', 'fig-l10.tex',
                'fig-l100.tex', 'fig-ms-ftab.tex', 'fig-lftab.tex')
# Locate groups of the freeze manifest, by cohort.
LOCATE_GROUPS = {'human_short': 'human_short-locate-6mb', 'human_pacbio': 'human_pacbio-locate-20mb',
                 'human_ont_modern': 'human_ont_modern-locate-10mb', 'guccione_skin': 'guccione_skin-locate-500kb'}
# Fig ops runs paper_figures.py reads, and the freeze group that holds each one's timed cells.
OPS_GROUPS = {OPS_RUNS['mtb152']: 'fig_ops-mtb152-heldout16-p100k-s600-v1',
              OPS_RUNS['human472']: 'fig_ops-human472-apr008-hap1-p100k-s600-v1',
              OPS2_RUNS['pacbio']['mtb152']: 'fig_ops-mtb152-pacbio-60mb-v2',
              OPS2_RUNS['pacbio']['human472']: 'fig_ops-human472-apr008-pacbio-60mb-v2'}
# Paper figure name: TikZ source in figures/tikz (as the paper's Makefile maps them). Fig rep is left
# out: its paper data come from runs outside the freeze manifests (results/README.md).
FIGMAP = {'ops': 'fig-ops-row', 'ms22_100': 'fig-ms22-100', 'ms22_10': 'fig-ms22-10', 'ms_oracle': 'fig-ms-oracle',
          'locate1': 'fig-l1', 'locate10': 'fig-l10', 'locate100': 'fig-l100', 'ms_ftab': 'fig-ms-ftab',
          'locate_ftab': 'fig-lftab'}


def run(*argv):
    print('+', ' '.join(map(str, argv)), flush=True)
    subprocess.run([sys.executable, *map(str, argv)], check=True, cwd=REPO)


def merge_ops(summaries, run_name, group):
    """A Fig ops run: the timed group's cells, each with ops_per_base from its -counts cell."""
    timed = json.loads((summaries/group/'summary.json').read_text())
    counts = {c['id']: c for c in json.loads((summaries/(group + '-counts')/'summary.json').read_text())['cells']}
    for c in timed['cells']:
        k = counts.get(c['id'])
        if k and k.get('status') == 'checked':
            o = k['ops_per_base']
            c['ops_per_base'] = o.get('ops_per_base', o)
        else:
            print('no counts for', group, c['id'])
    return timed


def stage(summaries, folder):
    for group in summaries.iterdir():
        if group.is_dir():
            (folder/group.name).symlink_to(group.resolve())
    for cohort, run_name in LOCATE_RUNS.items():
        (folder/run_name).symlink_to((summaries/LOCATE_GROUPS[cohort]).resolve())
    for run_name, group in OPS_GROUPS.items():
        (folder/run_name).mkdir()
        (folder/run_name/'summary.json').write_text(json.dumps(merge_ops(summaries, run_name, group), indent=1))


def compile_pdfs(out):
    tikz = out/'tikz'
    for name in TIKZ_SOURCES:
        if not (tikz/name).exists() or (tikz/name).resolve() != (SOURCES/name).resolve():
            shutil.copyfile(SOURCES/name, tikz/name)
    pdf = out/'pdf'
    pdf.mkdir(exist_ok=True)
    for paper, source in FIGMAP.items():
        subprocess.run(['latexmk', '-pdf', '-interaction=nonstopmode', '-outdir=build', f'{source}.tex'], cwd=tikz,
                       check=True, capture_output=True)
        shutil.copyfile(tikz/'build'/f'{source}.pdf', pdf/f'{paper}.pdf')
    print('wrote', pdf)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--results', default='results')
    ap.add_argument('--out', default='figures')
    ap.add_argument('--manifest', default='freeze-2885ee8-00e665de')
    ap.add_argument('--no-pdf', action='store_true')
    args = ap.parse_args()
    results, out = Path(args.results).resolve(), Path(args.out).resolve()
    summaries = out/'summaries'
    if summaries.exists():
        shutil.rmtree(summaries)
    run('-m', 'ms2exp.collect', '--manifest', results/args.manifest/'manifest.json', '--results-root', results,
        '--out', summaries, '--platform', 'aws', '--no-compare')
    # The run directories collect links for Fig rep point into results/; they are rebuilt on every run.
    rep_runs = summaries/'fig_rep'
    with tempfile.TemporaryDirectory() as tmp:
        stage(summaries, Path(tmp))
        run('-m', 'ms2exp.figures.paper_figures', '--ms-dir', tmp, '--chr22-dir', tmp, '--locate-dir', tmp,
            '--ops-dir', tmp, '--out', out/'draft')
    data = out/'tikz'/'data'
    run('-m', 'ms2exp.figures.paper_figures_tikz', '--csv', out/'draft/plotted-points.csv',
        '--captions', out/'draft/CAPTIONS.md', '--out', data)
    run('-m', 'ms2exp.figures.gap_annotations', data)
    run('-m', 'ms2exp.figures.fig_rep_tikz', '--a-dir', rep_runs, '--b-dir', rep_runs, '--out', data/'fig-rep',
        '--design', '4')
    shutil.rmtree(rep_runs)
    (out/'tex').mkdir(exist_ok=True)
    run('-m', 'ms2exp.figures.structure_bits', out/'tex/structure_bits.tex')
    if not args.no_pdf and shutil.which('latexmk'):
        compile_pdfs(out)


if __name__ == '__main__':
    main()
