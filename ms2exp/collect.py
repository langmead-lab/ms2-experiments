#!/usr/bin/env python3
"""Assemble freeze results into the per-group summaries ms2exp.figures.paper_figures reads, and check them.

Usage:
  python -m ms2exp.collect --manifest M --results-root DIR --out DIR [--platform aws|rockfish] [--no-compare]

For each cell group of the manifest, writes OUT/GROUP/summary.json in the shape of the
exploratory runs' summaries: top-level config.workload (reads, bases), config.orbit.revision,
config.rope.batch_size, and a cells list whose entries carry id, status, median_query_seconds,
median_load_seconds, peak_rss_kib, mbases_per_second, ns_per_base, locations_per_second,
layout, mode, interleave and ops_per_base. ms2exp.figures.paper_figures then reads OUT with its --ms-dir,
--chr22-dir, --locate-dir and --ops-dir pointed there and its run names mapped to the groups
(path changes only).

It also writes OUT/report.json and prints a report:
- cells without a done marker, and cells not checked (capped, timeout);
- locate cross-tool agreement: for every read set and cap, each ms2 cell's canonical output
  against each ropebwt3 cell's (ms2exp.validate.offsets.compare_files), which must agree. The
  canonical files stay in S3 rather than in results/, so --no-compare skips this step and the
  report keeps the comparisons of the results' own collect run;
- the D9 gate: MS dual on the full + seqids index against the plain full index, both on Rockfish,
  which must agree within NOISE of the plain cell's median query time;
- every result's ms2 commit equals the manifest's;
- the cross-thread check: cells that record output checksums (the thread-scaling matrix), grouped
  by tool, layout, mode, ropebwt3 length and read set, must have one output sha256 over every run
  and thread count.

A1 cells (kind rep) are also linked as OUT/fig_rep/runs/<name>, the run directories that
fig_rep_violin.py and fig_rep_tikz.py --design 2 read.
"""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

from ms2exp.lib import atomic_json, cell_key, file_sha256, load_manifest
from ms2exp.validate.offsets import compare_files

NOISE = 0.05
KEEP = ('id', 'status', 'layout', 'mode', 'interleave', 'block_bytes', 'cap', 'min_len', 'shift', 'batch',
        'median_query_seconds', 'median_load_seconds', 'query_seconds', 'load_seconds', 'peak_rss_kib',
        'mbases_per_second', 'ns_per_base', 'locations', 'smems', 'locations_per_second', 'ops_per_base',
        'query_spread_fraction', 'cap_s', 'tool', 'kind', 'figures', 'cache', 'host', 'platform', 'job', 'cell',
        'rep', 'threads', 'binding', 'output_sha256')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest', required=True)
    p.add_argument('--results-root', required=True, help='local copy of the results tree')
    p.add_argument('--out', required=True)
    p.add_argument('--platform', choices=('aws', 'rockfish'))
    p.add_argument('--no-compare', action='store_true', help='skip the locate cross-tool comparison')
    args = p.parse_args()
    manifest = load_manifest(args.manifest)
    base = Path(args.results_root)/manifest['id']
    out = Path(args.out)
    groups, missing, unchecked, wrong_commit, summaries = defaultdict(list), [], [], [], {}
    for job in manifest['jobs']:
        if args.platform and args.platform not in job['platforms']:
            continue
        for cell in job['cells']:
            if args.platform and args.platform not in cell.get('platforms', job['platforms']):
                continue
            folder = base/job['id']/'cells'/cell_key(cell['id'])
            done = folder.parent/(cell_key(cell['id']) + '.done')
            if not done.exists():
                missing.append(cell['id'])
                continue
            summary = json.loads((folder/'summary.json').read_text())
            if json.loads(done.read_text())['summary_sha256'] != file_sha256(folder/'summary.json'):
                raise SystemExit('Summary differs from its done marker: ' + cell['id'])
            if summary['ms2']['commit'] != manifest['ms2']['commit']:
                wrong_commit.append(cell['id'])
            if summary['status'] != 'checked':
                unchecked.append(f"{cell['id']}: {summary['status']}")
            summary['_folder'] = folder
            summaries[cell['id']] = summary
            groups[cell['group']].append(summary)
    for group, cells in sorted(groups.items()):
        first = cells[0]
        batches = sorted({c.get('batch') for c in cells if c['tool'] == 'ropebwt3' and c.get('batch')})
        doc = dict(format='ms2-freeze-group-v1', group=group, manifest_id=manifest['id'],
                   config=dict(workload=dict(reads=first['reads'], bases=first['bases'], reference=first['reference'],
                                             read_set=first['read_set']),
                               orbit=dict(revision=manifest['ms2']['commit'], repo=manifest['ms2']['repo']),
                               rope=dict(batch_size=batches[0] if len(batches) == 1 else None,
                                         binary_sha256=manifest['ropebwt3']['binary_sha256']),
                               protocol=manifest['protocol']),
                   cells=[{k: c[k] for k in KEEP if k in c} for c in sorted(cells, key=lambda c: c['id'])])
        atomic_json(out/group/'summary.json', doc)
    # A1: the run directories fig_rep_violin.py and fig_rep_tikz.py read, as OUT/fig_rep/runs/<name>.
    for cid, s in summaries.items():
        if s['kind'] == 'rep' and s['status'] == 'checked':
            link = out/'fig_rep'/'runs'/s['id']
            link.parent.mkdir(parents=True, exist_ok=True)
            if link.is_symlink():
                link.unlink()
            link.symlink_to((s['_folder']/s.get('rep_dir', 'rep')).resolve())
    # Locate: every ms2 cell against every ropebwt3 cell of the same read set and cap.
    locate = defaultdict(lambda: dict(ms2=[], rope=[]))
    for cid, s in summaries.items():
        if s['kind'] == 'locate' and s['status'] == 'checked':
            locate[(s['read_set'], s['cap'])]['ms2' if s['tool'] == 'ms2' else 'rope'].append(s)
    comparisons = []
    for (rset, cap), sides in sorted(locate.items()) if not args.no_compare else ():
        for a in sides['ms2']:
            for b in sides['rope']:
                ca, cb = a['_folder']/'canonical.gz', b['_folder']/'canonical.gz'
                result = compare_files((ca, a['check']['canonical']['sha256']),
                                       (cb, b['check']['canonical']['sha256']), cap, examples=3)
                comparisons.append(dict(read_set=rset, cap=cap, ms2=a['cell'], rope=b['cell'], status=result['status'],
                                        disagreements=result['disagreements'], examples=result['examples']))
    # D9 gate.
    d9 = []
    for cid, s in summaries.items():
        if cid.endswith('/dual-ordinary-plainfull') and s['status'] == 'checked':
            other = summaries.get(cid.rsplit('/', 1)[0] + '/dual-ordinary')
            if other and other['status'] == 'checked':
                ratio = other['median_query_seconds'] / s['median_query_seconds']
                d9.append(dict(group=s['group'], seqids_over_plain=ratio, passes=abs(ratio - 1) <= NOISE))
    threads = cross_thread(summaries.values())
    report = dict(manifest_id=manifest['id'], cells=len(summaries), missing=missing, unchecked=unchecked,
                  wrong_commit=wrong_commit, locate_comparisons='not run' if args.no_compare else comparisons,
                  d9_dual_seqids=d9, cross_thread=threads,
                  groups=sorted(groups))
    atomic_json(out/'report.json', report)
    print(f"{len(summaries)} cells with done markers in {len(groups)} groups; {len(missing)} without")
    for line in unchecked:
        print('  not checked:', line)
    for c in comparisons:
        print(f"  locate {c['read_set']} cap{c['cap']}: {c['ms2'].split('/')[-1]} vs {c['rope'].split('/')[-1]}: "
              f"{c['status']}")
    for d in d9:
        print(f"  D9 {d['group']}: seqids/plain = {d['seqids_over_plain']:.3f} ({'pass' if d['passes'] else 'FAIL'})")
    for t in threads:
        print(f"  threads {t['read_set']} {t['config']}: T = {t['threads']}: "
              f"{'one output' if t['passes'] else 'OUTPUTS DIFFER'}")
    bad = ([c for c in comparisons if c['status'] != 'agree'] + [d for d in d9 if not d['passes']] +
           [t for t in threads if not t['passes']])
    if wrong_commit:
        print('  results from another ms2 commit:', wrong_commit)
    return 1 if (bad or wrong_commit) else 0


def cross_thread(summaries):
    """One entry per tool configuration and read set whose checked cells record output checksums."""
    groups = defaultdict(list)
    for s in summaries:
        if s['status'] == 'checked' and s.get('output_sha256'):
            config = (s['tool'] + '-' + str(s['min_len']) if s['tool'] == 'ropebwt3'
                      else f"{s['tool']}-{s.get('layout')}-{s.get('mode')}")
            groups[(s['read_set'], config)].append(s)
    out = []
    for (rset, config), cells in sorted(groups.items()):
        digests = sorted({d for c in cells for d in c['output_sha256']})
        out.append(dict(read_set=rset, config=config, threads=sorted(c.get('threads', 1) for c in cells),
                        cells=sorted(c['cell'] for c in cells), sha256=digests, passes=len(digests) == 1))
    return out


if __name__ == '__main__':
    sys.exit(main())
