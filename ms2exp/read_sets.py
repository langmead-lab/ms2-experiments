#!/usr/bin/env python3
"""Records of the read sets a manifest queries: data/read-sets.json.

Usage:
  python -m ms2exp.read_sets --manifest results/freeze-2885ee8-00e665de/manifest.json
                             (--root DIR | --s3 s3://bucket/prefix [--profile P]) [--out data/read-sets.json]

For every reads/REFERENCE/SET/reads.fa in the manifest's inputs catalog, this reads the file (under
DIR, or streamed from the S3 prefix with `aws s3 cp URL -`, nothing written to disk), checks its
sha256 and size against the catalog, and records the number of reads, the bases and the read-length
distribution (minimum, median, N50, maximum). It adds what the files cannot show from SOURCES below:
where each set's reads come from and how simulated sets were made, as described in
ms2exp/build/d13_prep.py, ms2exp/sim/chr22_sim.py, ms2exp/figures/fig_rep.py and the paper.

The figures driver writes CAPTIONS.md from these records, and ms2exp.figures.paper_numbers checks
the paper's read-set numbers against them.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from ms2exp.lib import atomic_json

FORMAT = 'ms2-read-sets-v1'
HIFI = ('pbsim3 3.0.5 multi-pass subreads (QSHMM-RSII model, 10 passes) made into HiFi reads by ccs 6.4.0 '
        '(at least 3 passes, predicted accuracy at least 0.99), read lengths drawn from the real PacBio HiFi '
        'set, both strands')
CHR22 = 'APR008 haplotype 1 (Arab Pangenome Reference, GCA_050494745.1), long arm of chromosome 22'
APR008 = 'APR008 haplotype 1 (Arab Pangenome Reference, GCA_050494745.1), absent from HPRC2'
MTB_HELDOUT = 'the 8 mtb152 genomes held out of the index'
# Source and recipe of each set, keyed by the set's directory name with any -vN suffix removed; the
# first matching prefix wins.
SOURCES = [
    ('human_short-', dict(kind='real', platform='Illumina HiSeq 2500, 2x125 bp',
                          source='CHM1 hydatidiform mole cell line, SRR3099549; each end a separate query')),
    ('human_pacbio-', dict(kind='real', platform='PacBio Sequel II HiFi', source='HPRC sample HG00099, SRR26545347')),
    ('human_ont_modern-', dict(kind='real', platform='Oxford Nanopore R10.4.1',
                               source='HG002, ONT gm24385_2023.12 quality-passing reads')),
    ('guccione_skin-', dict(kind='real', platform='Illumina NovaSeq 6000, pairs of up to 2x151 bp',
                            source='atopic dermatitis skin swab, ERR14066769 (Guccione et al.)')),
    ('human472-apr008-hap1-p100k-s600', dict(kind='genome pieces', source=APR008,
                                            recipe='a seeded sample of 600 pieces of 100 kbases')),
    ('human472-apr008-pacbio-60mb', dict(kind='simulated', platform='PacBio HiFi', source=APR008, recipe=HIFI)),
    ('mtb152-heldout16-p100k-s600', dict(kind='genome pieces', source=MTB_HELDOUT + ' and their reverse complements',
                                         recipe='a seeded sample of 600 pieces of 100 kbases')),
    ('mtb152-pacbio-60mb', dict(kind='simulated', platform='PacBio HiFi', source=MTB_HELDOUT, recipe=HIFI)),
    ('mtb152-A2-exact', dict(kind='simulated', platform='PacBio-like', source=MTB_HELDOUT,
                             recipe='error-free templates with lengths drawn from the real PacBio HiFi set')),
    ('mtb152-A2-e', dict(kind='simulated', platform='PacBio-like', source=MTB_HELDOUT,
                         recipe='templates with lengths drawn from the real PacBio HiFi set and independent '
                                'substitutions, insertions and deletions at the rate in the set name')),
    ('simshort-60mb', dict(kind='simulated', platform='Illumina short reads', source=CHR22,
                           recipe='ART 2.5.8, HiSeq 2500 profile (-ss HS25), single-end 125 bases')),
    ('simpacbio-60mb', dict(kind='simulated', platform='PacBio HiFi', source=CHR22, recipe=HIFI)),
    ('pacbio-60mb', dict(kind='simulated', platform='PacBio HiFi', source=CHR22, recipe=HIFI)),
    ('simont-60mb', dict(kind='simulated', platform='Oxford Nanopore', source=CHR22,
                         recipe='pbsim3 3.0.5 standard ONT model (ERRHMM-ONT-HQ), 96% mean accuracy, read lengths '
                                'drawn from the real Nanopore set, both strands')),
]


def describe(set_name):
    stem = set_name.rsplit('-v', 1)[0]
    for prefix, record in SOURCES:
        if stem.startswith(prefix):
            return dict(record)
    raise KeyError(f'no source recorded for read set {set_name}; add it to SOURCES')


def fasta_stats(stream):
    """sha256, size and record lengths of a FASTA read from a binary stream."""
    digest, size, lengths, current = hashlib.sha256(), 0, [], None
    for line in stream:
        digest.update(line)
        size += len(line)
        if line.startswith(b'>'):
            if current is not None:
                lengths.append(current)
            current = 0
        else:
            current += len(line.rstrip(b'\r\n'))
    if current is not None:
        lengths.append(current)
    return digest.hexdigest(), size, lengths


def length_summary(lengths):
    s = sorted(lengths)
    half, run = sum(s) / 2, 0
    n50 = 0
    for length in reversed(s):
        run += length
        if run >= half:
            n50 = length
            break
    n = len(s)
    twice = 2 * s[n // 2] if n % 2 else s[n // 2 - 1] + s[n // 2]
    return dict(min=s[0], median=twice // 2 if twice % 2 == 0 else twice / 2, n50=n50, max=s[-1])


def read_one(name, args):
    if args.root:
        with open(Path(args.root)/name, 'rb') as fh:
            return fasta_stats(fh)
    cmd = ['aws', 's3', 'cp', f"{args.s3.rstrip('/')}/{name}", '-']
    if args.profile:
        cmd[1:1] = ['--profile', args.profile]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    out = fasta_stats(proc.stdout)
    if proc.wait():
        raise RuntimeError(f'{" ".join(cmd)} exited {proc.returncode}')
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--manifest', action='append', required=True)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument('--root', help='a directory holding the inputs under their logical names')
    src.add_argument('--s3', help='the inputs prefix, for example s3://ms2-freeze-159168350739/freeze/inputs')
    ap.add_argument('--profile')
    ap.add_argument('--out', default='data/read-sets.json')
    args = ap.parse_args()
    catalog = {}
    for path in args.manifest:
        for name, entry in json.loads(Path(path).read_text())['inputs'].items():
            if name.startswith('reads/') and name.endswith('/reads.fa'):
                catalog[name] = entry
    sets = {}
    for name in sorted(catalog):
        logical = name[len('reads/'):-len('/reads.fa')]
        sha, size, lengths = read_one(name, args)
        want = catalog[name]
        if (sha, size) != (want['sha256'], want['bytes']):
            sys.exit(f'{name}: sha256 {sha} and size {size} differ from the catalog ({want["sha256"]}, {want["bytes"]})')
        sets[logical] = dict(describe(logical.split('/')[1]), reads_fa=dict(sha256=sha, bytes=size),
                             reads=len(lengths), bases=sum(lengths), length=length_summary(lengths))
        print(f'{logical}: {len(lengths):,} reads, {sum(lengths):,} bases', flush=True)
    atomic_json(Path(args.out), dict(format=FORMAT, manifests=[json.loads(Path(p).read_text())['id']
                                                                for p in args.manifest], sets=sets))
    print('wrote', args.out)


if __name__ == '__main__':
    main()
