"""Checks of the committed input lists in data/lists against their sources.

human472-haplotypes.txt names the 472 AGC sets of human472.agc in archive order, with the PanSN
sample#haplotype prefix of each. check_haplotypes() compares that list with the sequence names of
Li's human472.fmd.len.gz, whose first two PanSN fields give one prefix per haplotype.
chr22-panel.tsv lists the chr22 components of the first 100 haplotypes of the donor order;
check_panel() checks its order against hprc2-donor-order.txt and the AGC set list.

Usage:
  python -m ms2exp.prep.lists haplotypes [--len human472.fmd.len.gz]
  python -m ms2exp.prep.lists panel
"""
import argparse
import gzip
from pathlib import Path
import sys

LISTS = Path(__file__).resolve().parents[2]/'data/lists'
SUFFIX = {'pri': '0', 'pat': '1', 'mat': '2', 'hap1': '1', 'hap2': '2'}


def pansn(agc_set):
    """PanSN prefix of an AGC set name, for example 200001_HG002.pat gives HG002#1."""
    sample, haplotype = agc_set.split('_', 1)[1].rsplit('.', 1)
    return f'{sample}#{SUFFIX[haplotype]}'


def table(path):
    """Rows of a tab-separated file with a header line, as dicts."""
    lines = Path(path).read_text().splitlines()
    header = lines[0].split('\t')
    return [dict(zip(header, line.split('\t'))) for line in lines[1:]]


def haplotypes(path=LISTS/'human472-haplotypes.txt'):
    rows = table(path)
    for row in rows:
        if pansn(row['agc_set']) != row['pansn_prefix']:
            raise ValueError('PanSN prefix differs from the set name: ' + row['agc_set'])
    names = [row['agc_set'] for row in rows]
    if len(set(names)) != len(names) or len({row['pansn_prefix'] for row in rows}) != len(rows):
        raise ValueError('Duplicate haplotype')
    return rows


def len_prefixes(path):
    """Distinct sample#haplotype prefixes of a ropebwt3 .len.gz file, in first-seen order."""
    seen = {}
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        for line in f:
            name = line.split('\t', 1)[0]
            fields = name.split('#')
            if len(fields) < 3:
                raise ValueError('Sequence name is not PanSN: ' + name)
            seen.setdefault('#'.join(fields[:2]), None)
    return list(seen)


def check_haplotypes(len_file, path=LISTS/'human472-haplotypes.txt', expected=472):
    rows = haplotypes(path)
    if expected is not None and len(rows) != expected:
        raise ValueError(f'{len(rows)} haplotypes, not {expected}')
    want = [row['pansn_prefix'] for row in rows]
    got = len_prefixes(len_file)
    if got != want:
        missing = sorted(set(want) - set(got))
        extra = sorted(set(got) - set(want))
        raise ValueError(f'Haplotypes differ from {len_file}: missing {missing[:5]}, extra {extra[:5]}, '
                         f'order equal {sorted(got) == sorted(want)}')
    return len(rows)


def donor_order(path=LISTS/'hprc2-donor-order.txt', count=100):
    """The first count haplotypes of the donor order, CHM13 and GRCh38 excluded (bench/ladder/plan.py)."""
    names = [n for n in Path(path).read_text().split() if 'CHM13' not in n and 'GRCh38' not in n]
    if len(set(names)) != len(names):
        raise ValueError('Duplicate donor in the order')
    return names[:count]


def check_panel(panel=LISTS/'chr22-panel.tsv', order=LISTS/'hprc2-donor-order.txt',
                sets=LISTS/'human472-haplotypes.txt', count=100):
    rows = table(panel)
    order = donor_order(order, count)
    known = {row['agc_set'] for row in haplotypes(sets)}
    if not set(order) <= known:
        raise ValueError('Donor order names a set absent from the archive')
    ranks = [int(row['rank']) for row in rows]
    if ranks != sorted(ranks):
        raise ValueError('Panel rows are not in donor order')
    for row in rows:
        r = int(row['rank'])
        if not 1 <= r <= count or order[r - 1] != row['haplotype']:
            raise ValueError(f"Rank {r} is not {row['haplotype']} in the donor order")
        if not row['contig'].startswith(pansn(row['haplotype']) + '#'):
            raise ValueError('Contig belongs to another haplotype: ' + row['contig'])
        if not 0 <= int(row['begin']) < int(row['end']) or row['orientation'] not in '+-':
            raise ValueError('Bad component: ' + row['contig'])
    if sorted(set(ranks)) != list(range(1, count + 1)):
        raise ValueError('A haplotype of the panel has no chr22 component')
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='what', required=True)
    h = sub.add_parser('haplotypes')
    h.add_argument('--len', help='Li\'s human472.fmd.len.gz; without it only the list itself is checked')
    sub.add_parser('panel')
    args = p.parse_args(argv)
    if args.what == 'haplotypes':
        n = check_haplotypes(args.len) if args.len else len(haplotypes())
        print(f'{n} haplotypes checked')
    else:
        rows = check_panel()
        print(f'{len(rows)} components of {len({r["haplotype"] for r in rows})} haplotypes checked')


if __name__ == '__main__':
    sys.exit(main())
