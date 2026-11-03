#!/usr/bin/env python3
"""The chr22 references behind every chr22 index: chr22 of the first N haplotypes of HPRC2, from human472.agc.

data/lists/chr22-panel.tsv lists, for the first 100 haplotypes of data/lists/hprc2-donor-order.txt (CHM13 and
GRCh38 excluded), each contig interval that RagTag v2.1.0 (ragtag.py scaffold against CHM13, the lab's runs
in data/evidence/ragtag/) placed on CHM13 chr22. The table is the chr22 part of the frozen ladder plan of
ms-bakery bench/ladder/plan.py; the chr22-hN reference is the rows of rank 1 to N.

For one N, this extracts each component with agc getctg, uppercases it, splits it at every non-ACGT run into
records named h<haplotype number>_c<component>_s<segment>, and writes in OUT/chr22-hN:
  reference.fa         the forward segments (the ropebwt3 index input; ropebwt3 build -d adds the other strand);
  reference-both.txt   the forward segments then their reverse complements, joined by N (MONI, PHONI, AugPHONI);
  reference-both.fa    the same text as one FASTA record;
  manifest.json        components, extracted-FASTA sha256s, segment coordinates and reference.fa's sha256.
Components stay separate; no RagTag gaps or joins enter the text. --check compares reference.fa's sha256 with
the frozen manifest in data/evidence/rockfish/data/ladder/chr22-hN/manifest.json.

--from-agp DIR rebuilds the panel rows from RagTag AGP files (DIR/<haplotype>/ragtag_output/ragtag.scaffold.agp,
with haplotype directories named as in the lab's runs, for example HG002.1) and reports any difference from
chr22-panel.tsv. The AGPs are not public, so this is a check for the lab, not a step of the rebuild.

Usage:
  python3 data/prep/chr22_panel.py extract N --agc human472.agc --agc-binary agc --out DIR [--check]
  python3 data/prep/chr22_panel.py from-agp --agp-root DIR
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from ms2exp.prep import lists  # noqa: E402

PANEL = REPO/'data/lists/chr22-panel.tsv'
EVIDENCE = REPO/'data/evidence/rockfish/data/ladder'
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def records(path):
    name, parts = None, []
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(parts)
                name, parts = line[1:].strip(), []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, ''.join(parts)


def components(count, panel=PANEL):
    rows = lists.check_panel(panel)
    return [dict(haplotype=r['haplotype'], contig=r['contig'], begin=int(r['begin']), end=int(r['end']),
                 orientation=r['orientation'], contig_md5=r['contig_md5'])
            for r in rows if int(r['rank']) <= count]


def extract(count, archive, agc, out, cache=None):
    out = Path(out)/f'chr22-h{count}'
    out.mkdir(parents=True, exist_ok=False)
    cache = Path(cache or out/'components')
    cache.mkdir(parents=True, exist_ok=True)
    metadata, total = [], 0
    with (out/'reference.fa').open('w') as dst:
        for i, r in enumerate(components(count)):
            key = hashlib.sha256(json.dumps(r, sort_keys=True).encode()).hexdigest()
            path = cache/(key + '.fa')
            if not path.exists():
                tmp = path.with_suffix('.tmp')
                region = f"{r['contig']}@{r['haplotype']}:{r['begin']}-{r['end'] - 1}"
                with tmp.open('w') as f:
                    subprocess.run([agc, 'getctg', '-t', '2', archive, region], stdout=f, check=True)
                seq = ''.join(s for _, s in records(tmp)).upper()
                if len(seq) != r['end'] - r['begin']:
                    raise ValueError('AGC length mismatch: ' + str(r))
                tmp.replace(path)
            seq = ''.join(s for _, s in records(path)).upper()
            segments = []
            for j, m in enumerate(re.finditer('[ACGT]+', seq)):
                name = f"h{r['haplotype'].split('_')[0]}_c{i}_s{j}"
                dst.write(f'>{name}\n{m.group()}\n')
                total += len(m.group())
                segments.append(dict(record=name, begin=r['begin'] + m.start(), end=r['begin'] + m.end()))
            metadata.append(dict(r, extracted_fasta_sha256=digest(path), segments=segments))
    first = True
    with (out/'reference-both.txt').open('w') as txt, (out/'reference-both.fa').open('w') as fa:
        fa.write('>reference\n')
        for reverse in (False, True):
            for _, seq in records(out/'reference.fa'):
                if reverse:
                    seq = seq.translate(COMPLEMENT)[::-1]
                sep = '' if first else 'N'
                txt.write(sep + seq)
                fa.write(sep + seq)
                first = False
        fa.write('\n')
    manifest = dict(id=f'chr22-h{count}', haplotypes=count, archive=str(archive), panel=str(PANEL.relative_to(REPO)),
                    panel_sha256=digest(PANEL), sources=metadata, forward_acgt_bases=total,
                    reference_sha256=digest(out/'reference.fa'), both_strands=True,
                    boundary_policy='Native mapped component intervals remain separate; split non-ACGT runs. '
                                    'No synthetic RagTag joins. Search both strands.')
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def frozen_sha256(count):
    path = EVIDENCE/f'chr22-h{count}'/'manifest.json'
    return json.loads(path.read_text())['reference_sha256'] if path.exists() else None


def agp_rows(path, chromosome='chr22'):
    """W lines of a RagTag AGP whose object is the CHM13 chromosome, as (contig, begin, end, orientation)."""
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        f = line.split('\t')
        m = re.search(r'#(chr(?:[0-9]+|X|Y|M))_RagTag$', f[0])
        if f[4] == 'W' and m and m[1] == chromosome:
            rows.append((f[5], int(f[6]) - 1, int(f[7]), f[8]))
    return rows


def agp_directory(haplotype):
    """The lab's RagTag run directory name of an AGC set, for example 200001_HG002.pat gives HG002.1."""
    return lists.pansn(haplotype).replace('#', '.')


def from_agp(root, panel=PANEL):
    rows = lists.check_panel(panel)
    differences = []
    for h in dict.fromkeys(r['haplotype'] for r in rows):
        agp = Path(root)/agp_directory(h)/'ragtag_output/ragtag.scaffold.agp'
        want = [(r['contig'], int(r['begin']), int(r['end']), r['orientation']) for r in rows if r['haplotype'] == h]
        got = agp_rows(agp)
        if got != want:
            differences.append(dict(haplotype=h, agp=str(agp), panel=want, agp_rows=got))
        elif digest(agp) != next(r['agp_sha256'] for r in rows if r['haplotype'] == h):
            differences.append(dict(haplotype=h, agp=str(agp), note='same rows, different file sha256'))
    return differences


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='action', required=True)
    e = sub.add_parser('extract')
    e.add_argument('count', type=int)
    e.add_argument('--agc', required=True, help='human472.agc')
    e.add_argument('--agc-binary', default='agc')
    e.add_argument('--out', required=True)
    e.add_argument('--cache', help='shared directory of extracted components')
    e.add_argument('--check', action='store_true')
    a = sub.add_parser('from-agp')
    a.add_argument('--agp-root', required=True)
    args = p.parse_args(argv)
    if args.action == 'extract':
        m = extract(args.count, args.agc, args.agc_binary, args.out, args.cache)
        print(m['id'], m['forward_acgt_bases'], 'forward ACGT bases', m['reference_sha256'])
        if args.check:
            want = frozen_sha256(args.count)
            if want is None:
                print('No frozen manifest for this N')
            elif want != m['reference_sha256']:
                sys.exit(f'reference.fa differs from the frozen reference ({want})')
            else:
                print('reference.fa equals the frozen reference')
    else:
        d = from_agp(args.agp_root)
        print(json.dumps(d, indent=1) if d else 'Every haplotype matches its AGP')
        sys.exit(1 if d else 0)


if __name__ == '__main__':
    main()
