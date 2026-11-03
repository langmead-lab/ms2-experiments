#!/usr/bin/env python3
"""Compare two toy runs (group summaries, per-cell records and done markers), or one toy run with
the expected summaries in toy/expected.

Usage:
  compare.py TOY_A TOY_B [--collected-a NAME] [--collected-b NAME] [--results-a NAME] [--results-b NAME]
  compare.py --expected toy/expected TOY [--collected NAME] [--own-tools]

Fields that record time, memory, cache state or host are left out, as are paths under the toy
directories, staging directories and manifest ids. Everything else must be equal: statuses,
checks, operation counts, locations, commands and input sha256s. Index inputs built separately for
each toy differ in sha256 only through their header's build path and time, so cell records that
name them are reported apart from other differences. Exits 1 on any other difference.

With --expected, the reference is a directory of group summaries, report.json and
manifest-toy.json, as toy/expected holds, so that a run anywhere can be checked against the lab's.
The ms2 source location (config.orbit.repo and the manifest's ms2.repo) is left out too. The
manifests' protocol, seed, jobs, ms2 commit and ropebwt3 pin must be equal, and so must the sha256
of every input except those each toy builds with its own compiler or records itself: ms2/*,
tools/ms2-*, tools/ms_check and traces/*. That checks the fetched inputs and the files the toy
derives from them (the two cut read sets, the .len.gz and the SSA) byte for byte. --own-tools also
leaves out the ropebwt3 binary's sha256 and every tools/* input, for a run with ropebwt3 and the
legacy tools built from source.
"""
import argparse
import json
from pathlib import Path
import re
import sys

VOLATILE = re.compile(r'^(wall_seconds|native_wall_seconds|user_seconds|system_seconds|peak_rss_kib|cache|host|'
                      r'started|finished|created|at|load_seconds|query_seconds|median_query_seconds|'
                      r'median_load_seconds|query_spread_fraction|mbases_per_second|ns_per_base|'
                      r'locations_per_second|index_load_s|ftab_s|query_s|read_s|parse_s|format_s|total_s|'
                      r'last_loaded_s|last_processed_s|manifest_id|seconds|heartbeat)$')


def norm(x, roots):
    if isinstance(x, dict):
        return {k: norm(v, roots) for k, v in x.items() if not VOLATILE.match(k)}
    if isinstance(x, list):
        return [norm(v, roots) for v in x]
    if isinstance(x, str):
        for root in roots:
            x = x.replace(root, 'TOY')
        x = re.sub(r'/tmp/[^/]+', 'STAGE', x)
        x = re.sub(r'toy-[0-9a-f]{7}-[0-9a-f]{8}', 'MANIFEST', x)
    return x


def diff(a, b, path=''):
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b)):
            out += [f'{path}/{k}: missing on one side'] if k not in a or k not in b else diff(a[k], b[k], f'{path}/{k}')
        return out
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        return [d for i, (x, y) in enumerate(zip(a, b)) for d in diff(x, y, f'{path}[{i}]')]
    return [] if a == b else [f'{path}: {str(a)[:80]!r} != {str(b)[:80]!r}']


def load(path):
    return json.loads(Path(path).read_text())


# Inputs a toy builds or records itself; their bytes differ from one build to the next.
OWN_BUILD = ('ms2/', 'tools/ms2-', 'tools/ms_check', 'traces/')


def drop(x, keys):
    """X without the dict entries at the given (dotted) key paths."""
    x = json.loads(json.dumps(x))
    for path in keys:
        node = x
        *parents, leaf = path.split('.')
        for k in parents:
            node = node.get(k, {}) if isinstance(node, dict) else {}
        if isinstance(node, dict):
            node.pop(leaf, None)
    return x


def against_expected(expected, toy, collected, own_tools):
    roots = [str(expected), str(toy)]
    bad = []
    me, mt = load(expected/'manifest-toy.json'), load(toy/'manifest-toy.json')
    keys = ['ms2.repo'] + (['ropebwt3.binary_sha256'] if own_tools else [])
    for k in ('matrix', 'seed', 'protocol', 'ms2', 'ropebwt3', 'jobs'):
        d = diff(norm(drop({k: me.get(k)}, keys), roots), norm(drop({k: mt.get(k)}, keys), roots))
        bad += [f'manifest{line}' for line in d]
    skip = OWN_BUILD + (('tools/',) if own_tools else ())
    ie, it = me['inputs'], mt['inputs']
    names = sorted(n for n in set(ie) | set(it) if not n.startswith(skip))
    differ = [n for n in names if ie.get(n, {}).get('sha256') != it.get(n, {}).get('sha256')]
    print(f'manifest: {len(names)} inputs compared by sha256, {len(differ)} differ')
    bad += [f'input {n}: {ie.get(n, {}).get("sha256", "missing")[:16]} != {it.get(n, {}).get("sha256", "missing")[:16]}'
            for n in differ]
    groups = sorted(p.parent.name for p in expected.glob('*/summary.json'))
    if groups != sorted(p.parent.name for p in (toy/collected).glob('*/summary.json')):
        bad.append('group sets differ')
    keys = ['config.orbit.repo'] + (['config.rope.binary_sha256'] if own_tools else [])
    for g in groups:
        got = toy/collected/g/'summary.json'
        if not got.exists():
            continue
        d = diff(norm(drop(load(expected/g/'summary.json'), keys), roots), norm(drop(load(got), keys), roots))
        print(f'group {g}: {len(d)} differences')
        bad += [f'{g}{line}' for line in d]
    d = diff(norm(load(expected/'report.json'), roots), norm(load(toy/collected/'report.json'), roots))
    print('report.json:', d or 'equal')
    bad += d
    for line in bad[:40]:
        print('  ', line)
    print('equivalent' if not bad else f'{len(bad)} differences')
    return 1 if bad else 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('a', nargs='?', help='first toy directory (with --expected: the toy directory to check)')
    p.add_argument('b', nargs='?')
    p.add_argument('--expected', help='a directory of expected summaries, such as toy/expected')
    p.add_argument('--collected', default='collected', help='with --expected: the toy\'s collect directory name')
    p.add_argument('--own-tools', action='store_true', help='with --expected: ropebwt3 and legacy tools built here')
    p.add_argument('--collected-a', default='collected')
    p.add_argument('--collected-b', default='collected')
    p.add_argument('--results-a', default='results')
    p.add_argument('--results-b', default='results')
    args = p.parse_args()
    if args.expected:
        if not args.a or args.b:
            p.error('--expected EXPECTED takes one toy directory')
        return against_expected(Path(args.expected).resolve(), Path(args.a).resolve(), args.collected, args.own_tools)
    if not args.a or not args.b:
        p.error('give two toy directories, or --expected EXPECTED and one')
    a, b = Path(args.a).resolve(), Path(args.b).resolve()
    roots = [str(a), str(b)]
    ma, mb = load(a/'manifest-toy.json'), load(b/'manifest-toy.json')
    ia, ib = load(a/'toy-inputs.json')['inputs'], load(b/'toy-inputs.json')['inputs']
    rebuilt = {k for k in set(ia) | set(ib) if ia.get(k, {}).get('sha256') != ib.get(k, {}).get('sha256')}
    hashes = {ia[k]['sha256'] for k in rebuilt if k in ia} | {ib[k]['sha256'] for k in rebuilt if k in ib}
    print('inputs with different sha256:', sorted(rebuilt))
    skip = ('created', 'id', 'generator', 'inputs')
    bad = diff({k: v for k, v in ma.items() if k not in skip}, {k: v for k, v in mb.items() if k not in skip})
    print('manifest differences besides created, id, generator and inputs:', bad)
    ca, cb = a/args.collected_a, b/args.collected_b
    groups = sorted(p.parent.name for p in ca.glob('*/summary.json'))
    if groups != sorted(p.parent.name for p in cb.glob('*/summary.json')):
        bad.append('group sets differ')
    for g in groups:
        d = diff(norm(load(ca/g/'summary.json'), roots), norm(load(cb/g/'summary.json'), roots))
        print(f'group {g}: {len(d)} differences')
        bad += d
    d = diff(norm(load(ca/'report.json'), roots), norm(load(cb/'report.json'), roots))
    print('report.json:', d or 'equal')
    bad += d
    ra, rb = a/args.results_a/ma['id'], b/args.results_b/mb['id']
    files = sorted(p.relative_to(ra) for p in ra.glob('*/cells/*/*.json'))
    if files != sorted(p.relative_to(rb) for p in rb.glob('*/cells/*/*.json')):
        bad.append('cell record sets differ')
    named = 0
    for f in files:
        d = diff(norm(load(ra/f), roots), norm(load(rb/f), roots))
        if d and all(any(h in line for h in hashes) for line in d):
            named += 1
        elif d:
            bad += [f'{f}{line}' for line in d]
    print(f'{len(files)} cell records; {named} differ only in the sha256 of rebuilt inputs')
    done = sorted(p.relative_to(ra) for p in ra.glob('*/cells/*.done'))
    print(f'{len(done)} done markers')
    if done != sorted(p.relative_to(rb) for p in rb.glob('*/cells/*.done')):
        bad.append('done marker sets differ')
    for line in bad[:40]:
        print('  ', line)
    print('equivalent' if not bad else f'{len(bad)} differences')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
