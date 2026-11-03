#!/usr/bin/env python3
"""Seeded read pools of the real read sets: the first step behind the CHM1, PacBio and Nanopore sets.

This is the preparation that made the pools under real-hprc2-1d580b6/inputs/COHORT on Rockfish
(data/evidence/rockfish/results/real-hprc2-1d580b6/inputs/COHORT/complete.json records each run). The
frozen read sets were cut from those pools, which are kept as they are; this module is their provenance
and is not rerun by data/prep/run_d13.sh. Rerunning it with SeqKit 2.10.1 and the same sources should
give the same files, which complete.json lets one check.

For one cohort:
  1. seqkit sample -j 4 -s 20261004 -p FRACTION on each source FASTQ, FRACTION = min(1, 1e9 / source
     bases) for about 1 Gbase (1 for human_ont_modern, whose source is already the 0.8% samtools sample of
     data/download/prepare_hg002.py). Both mate files of a pair use one seed and fraction; the script
     checks that the selected mates come in the same order.
  2. Every selected read, uppercased, becomes record r0, r1, ... of selected.fa, mate 1 file first;
     original-read-names.tsv gives each record's mate, original name and length.
  3. seqkit sample -s 20261005 on selected.fa draws the nested sets core (p 0.2, about 200 Mbases),
     pilot and diagnostic (0.001) and offset (0.0001). Records are drawn independently, so the two mates
     of a pair are not kept together in these sets.
  4. Each set is split at every non-ACGT run into records q0, q1, ... (SET.fa), with SET.map.tsv (q name,
     selected record, start, end), SET.offsets.u64 and SET.one.fa; selected.fa itself is the scaling set.

Usage: python -m ms2exp.prep.pool COHORT --sources DIR --out DIR --seqkit PATH
  DIR for --sources holds the downloads in the layout of data/download (SRR3099549/..., SRR26545347/...,
  HG002-Q28-s008-seed20261003/sample.fastq.gz); --out gets COHORT/.
"""
import argparse
import array
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import re
import subprocess
import sys
import time

SEED = 20261004
NESTED_SEED = 20261005
TARGET = 1e9
# Source files and their published base counts (the denominator of the sampling fraction).
SETS = {
    'human_short': (['SRR3099549/SRR3099549_1.fastq.gz', 'SRR3099549/SRR3099549_2.fastq.gz'], 57581118000),
    'human_pacbio': (['SRR26545347/SRR26545347_subreads.fastq.gz'], 166559881890),
    'human_ont_modern': (['HG002-Q28-s008-seed20261003/sample.fastq.gz'], 1022850476),
}
NESTED = (('core', .2), ('pilot', .001), ('diagnostic', .001), ('offset', .0001))


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fastq(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)) as f:
        while header := f.readline():
            seq = f.readline().rstrip()
            plus = f.readline()
            qual = f.readline().rstrip()
            if not header.startswith('@') or not plus.startswith('+') or len(seq) != len(qual):
                raise ValueError('Malformed FASTQ record in ' + str(path))
            yield header[1:].split()[0], seq, qual


def fasta_records(path):
    """(name, sequence) of a FASTA file, sequences joined across lines (SeqKit wraps its FASTA output)."""
    with open(path) as f:
        name, parts = None, []
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
        if name is not None:
            yield name, ''.join(parts)


def mate_key(name):
    return re.sub(r'[/\.][12]$', '', name)


def check_mates(first, second):
    """The two mate samples hold the same pairs in the same order."""
    for a, b in itertools.zip_longest(fastq(first), fastq(second)):
        if a is None or b is None or mate_key(a[0]) != mate_key(b[0]):
            raise ValueError(f'Paired names differ: {a and a[0]} {b and b[0]}')


def write_selected(samples, d):
    """selected.fa and original-read-names.tsv from the sampled FASTQ files; returns the read count."""
    count = 0
    with (d/'selected.fa').open('w') as out, (d/'original-read-names.tsv').open('w') as ids:
        for mate, sample in enumerate(samples, 1):
            for name, seq, _ in fastq(sample):
                key = f'r{count}'
                count += 1
                out.write(f'>{key}\n{seq.upper()}\n')
                ids.write(f'{key}\t{mate}\t{name}\t{len(seq)}\n')
    return count


def split(source, d, label):
    """SET.fa, SET.map.tsv, SET.offsets.u64 and SET.one.fa: source records split at non-ACGT runs."""
    n = bases = ambiguous = affected = 0
    offsets = [0]
    out, mapping = d/(label + '.fa'), d/(label + '.map.tsv')
    with out.open('w') as dest, mapping.open('w') as mp:
        for key, s in fasta_records(source):
            runs = list(re.finditer('[ACGT]+', s))
            removed = len(s) - sum(len(x[0]) for x in runs)
            ambiguous += removed
            affected += bool(removed)
            for run in runs:
                dest.write(f'>q{n}\n{run[0]}\n')
                mp.write(f'q{n}\t{key}\t{run.start()}\t{run.end()}\n')
                n += 1
                bases += len(run[0])
                offsets.append(bases)
    if n < 2:
        raise ValueError(f'{label}: fewer than two records')
    a = array.array('Q', offsets)
    if sys.byteorder != 'little':
        a.byteswap()
    (d/(label + '.offsets.u64')).write_bytes(a.tobytes())
    with out.open() as inp:
        (d/(label + '.one.fa')).write_text(next(inp) + next(inp))
    return dict(reads=n, bases=bases, one_bases=offsets[1], ambiguous_bases_removed=ambiguous,
                affected_original_reads=affected, sha256=sha256(out), map_sha256=sha256(mapping))


def prepare(cohort, sources, out, seqkit):
    d = Path(out)/cohort
    d.mkdir(parents=True, exist_ok=True)
    if (d/'complete.json').exists():
        return json.loads((d/'complete.json').read_text())
    paths, total = SETS[cohort]
    prob = 1 if cohort == 'human_ont_modern' else min(1, TARGET/total)
    state = dict(cohort=cohort, status='running', started=time.time(), seed=SEED,
                 source_total_bases_for_fraction=total, source_fraction=prob,
                 seqkit_version=subprocess.check_output([str(seqkit), 'version'], text=True),
                 seqkit_sha256=sha256(seqkit), commands=[], sources=[])
    save(d/'state.json', state)
    try:
        samples = []
        for i, rel in enumerate(paths):
            source = Path(sources)/rel
            state['sources'].append(dict(path=str(source), bytes=source.stat().st_size, mtime=source.stat().st_mtime))
            dest = d/f'selected-{i}.fq.gz'
            argv = [str(seqkit), 'sample', '-j', '4', '-s', str(SEED), '-p', str(prob), str(source), '-o', str(dest)]
            state['commands'].append(argv)
            save(d/'state.json', state)
            if not (dest.exists() and (d/'selected.fa').exists()):
                subprocess.run(argv, check=True)
            samples.append(dest)
        if len(samples) == 2:
            check_mates(*samples)
        count = write_selected(samples, d)
        files = {'scaling': d/'selected.fa'}
        for label, p in NESTED:
            dest = d/(label + '.raw.fa')
            argv = [str(seqkit), 'sample', '-j', '4', '-s', str(NESTED_SEED), '-p', str(p), str(d/'selected.fa'),
                    '-o', str(dest)]
            state['commands'].append(argv)
            if not dest.exists():
                subprocess.run(argv, check=True)
            files[label] = dest
        stats = {label: split(f, d, label) for label, f in files.items()}
        state.update(status='complete', finished=time.time(), sets=stats, selected_original_reads=count)
        save(d/'complete.json', state)
        save(d/'state.json', state)
        return state
    except Exception as e:
        state.update(status='failed', error=repr(e))
        save(d/'state.json', state)
        raise


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cohort', choices=sorted(SETS))
    p.add_argument('--sources', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--seqkit', required=True, help='SeqKit 2.10.1 (sha256 ea5aef8d... in tools.lock.json prep_tools)')
    args = p.parse_args(argv)
    state = prepare(args.cohort, args.sources, args.out, Path(args.seqkit))
    print(json.dumps({k: v for k, v in state.items() if k in ('cohort', 'status', 'sets')}, indent=1))


if __name__ == '__main__':
    main()
