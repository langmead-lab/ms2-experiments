#!/usr/bin/env python3
"""The untimed screen that chose the skin metagenome run (ERR14066769) among the study's skin runs.

It ranks every deposited skin run of Guccione et al. (PRJEB83637) by residual human sequence: for a
seeded bottom-k sample of 10,000 pairs per run (both mates as independent reads), the fraction of
reads whose longest exact human472 match, from ropebwt3 mem -l1, is at least 50 bases. The run with
the highest fraction among runs with at least 60 Mbases of usable bases is chosen. The frozen result
is data/evidence/skin/guccione-screen-v1.json (77 runs, selected ERR14066769, no ties).

The skin runs come from data/download/download_skin.py --all, whose skin-verified.json lists each
file with its byte count, MD5 and sha256.

Usage: python -m ms2exp.prep.skin_screen --root DIR --ropebwt3 PATH --fmd human472.fmd --output JSON --work DIR
"""
import argparse
import gzip
import hashlib
import heapq
import json
import multiprocessing
import os
from fractions import Fraction
from pathlib import Path
import re
import socket
import subprocess
import time

from ms2exp.sim.chr22_sim import sha256

ROPEBWT3_SHA256 = 'cfad8475911f9298bf53177bf4a578a6d133f0464f4b9970bef9afb04697d53d'
INDEX_SHA256 = '10e89eb17f17484e063746b7515a6098ab1f3f398f19db94fed4e710e05c7626'
ACGT = re.compile('[ACGT]+')


def segments(sequence):
    """Uppercase, then split at every non-ACGT base, as the historical cohorts do."""
    return [(m.start(), m.end(), m.group()) for m in ACGT.finditer(sequence.upper())]


def usable(sequence):
    upper = sequence.upper()
    return sum(map(upper.count, 'ACGT'))


def mate_key(name):
    return name[:-2] if name.endswith(('/1', '/2')) else name


def runs(manifest, root, cohort='skin'):
    """Map each run to its verified R1 and R2 entries from the download manifest."""
    entries = [e for e in json.loads(Path(manifest).read_text()) if e['cohort'] == cohort]
    found = {}
    for e in entries:
        mate = {'_1.fastq.gz': 0, '_2.fastq.gz': 1}.get(e['relative_path'][-11:])
        if mate is None or Path(e['relative_path']).parts[:2] != (cohort, e['run']):
            raise ValueError('Unexpected deposited file: '+e['relative_path'])
        pair = found.setdefault(e['run'], [None, None])
        if pair[mate] is not None:
            raise ValueError('Duplicate deposited mate: '+e['relative_path'])
        pair[mate] = dict(path=str(Path(root)/e['relative_path']), bytes=e['bytes'], sha256=e['sha256'], md5=e['md5'])
    if any(None in pair for pair in found.values()):
        raise ValueError('Run lacks one mate file')
    return {run: tuple(found[run]) for run in sorted(found)}


class _Hashing:
    def __init__(self, path):
        self.handle = open(path, 'rb');self.digest = hashlib.sha256();self.bytes = 0

    def read(self, size=-1):
        block = self.handle.read(size)
        self.digest.update(block);self.bytes += len(block)
        return block

    def drain(self):
        while self.read(16 << 20):pass
        self.handle.close()


class Pairs:
    """Iterate (ordinal, name, R1, R2); finish() reads the rest and checks both files."""
    def __init__(self, r1, r2):
        self.entries = (r1, r2)
        self.raw = [_Hashing(e['path']) for e in self.entries]
        self.text = [gzip.open(raw, 'rt') for raw in self.raw]

    def _records(self, stream):
        while header := stream.readline():
            sequence = stream.readline().rstrip('\n');plus = stream.readline();quality = stream.readline().rstrip('\n')
            if not header.startswith('@') or not plus.startswith('+') or len(sequence) != len(quality):
                raise ValueError('Malformed FASTQ record')
            yield header[1:].split()[0], sequence

    def __iter__(self):
        ordinal = 0
        one, two = (self._records(t) for t in self.text)
        for a in one:
            b = next(two, None)
            if b is None or mate_key(a[0]) != mate_key(b[0]):
                raise ValueError('Mate files differ at pair '+str(ordinal))
            yield ordinal, mate_key(a[0]), a[1], b[1]
            ordinal += 1
        if next(two, None) is not None:
            raise ValueError('Mate files differ in length')

    def finish(self):
        for raw, entry in zip(self.raw, self.entries):
            raw.drain()
            if raw.bytes != entry['bytes'] or raw.digest.hexdigest() != entry['sha256']:
                raise ValueError('Deposited file differs from verified manifest: '+entry['path'])
        return [dict(path=e['path'], bytes=e['bytes'], sha256=e['sha256']) for e in self.entries]


SEED = 'guccione-skin-human472-screen-20261016'
PAIRS = 10000
THRESHOLD = 50
GRID = (19, 25, 31, 40, 50, 75, 100, 150)
MIN_BASES = 60000000
CAVEAT = ('Measures residual human sequence in deposited, host-filtered reads: the fraction of sampled reads whose '
        'longest exact Human472 match reaches the threshold. It is not the paper\'s pre-filtration host fraction '
        'and supports no host-fraction claim. Untimed; no timing from this screen is publishable.')
JUSTIFICATION = ('Random exact matches against both strands of 472 haplotypes (about 3e12 bases) have expected longest '
                 'length near log4(3e12 x 150), about 24 bases, for a 150-base read, and reach 31 by chance in roughly '
                 '1e-4 of reads. 50 bases is twice that expectation, so chance hits and short conserved or low-complexity '
                 'stretches rarely count, while a 150-base human read at Illumina error rates almost always keeps an '
                 'error-free stretch of 50 or more. Fractions at every grid threshold, including 31, are recorded.')


def key(seed, run, name):
    return int.from_bytes(hashlib.sha256(f'{seed}:{run}:{name}'.encode()).digest()[:8], 'big')


def sample_run(job):
    """Bottom-k pairs by seeded name hash: exact size, one pass, independent of file order."""
    run, r1, r2, pairs, seed = job
    heap = [];totals = dict(pairs=0, reads=0, bases=0, usable_bases=0, reads_without_acgt=0)
    stream = Pairs(r1, r2)
    for ordinal, name, a, b in stream:
        ua, ub = usable(a), usable(b)
        totals['pairs'] += 1;totals['reads'] += 2;totals['bases'] += len(a)+len(b)
        totals['usable_bases'] += ua+ub;totals['reads_without_acgt'] += (ua == 0)+(ub == 0)
        item = (-key(seed, run, name), -ordinal, name, a, b)
        if len(heap) < pairs:heapq.heappush(heap, item)
        elif item > heap[0]:heapq.heapreplace(heap, item)
    files = stream.finish()
    sample = sorted((-o, name, a, b) for _, o, name, a, b in heap)
    return dict(run=run, files=files, sampled_pairs=len(sample), **totals), sample


def default_runner(argv, stdout, stderr):
    with open(stdout, 'w') as out, open(stderr, 'w') as err:
        subprocess.run(argv, stdout=out, stderr=err, check=True,
                       env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))


def longest_matches(path, count):
    """Longest SMEM per query segment from ropebwt3 -l1 output; every segment must match."""
    best = [0]*count
    with open(path) as f:
        for line in f:
            fields = line.split()
            if not fields:continue
            name, start, end, hits = fields[:4]
            i = int(name[1:])
            if name != f'q{i}' or not 0 <= i < count:raise ValueError('Unexpected screen record '+name)
            start, end = int(start), int(end)
            if int(hits) > 0 and end-start > best[i]:best[i] = end-start
    if 0 in best:raise ValueError('Screen output lacks a match for some segment')
    return best


def summarize(values, threshold, grid=GRID):
    n = len(values);ordered = sorted(values)
    nearest = lambda q: ordered[min(n-1, max(0, -(-n*q//100)-1))]
    histogram = {}
    for v in ordered:histogram[str(v)] = histogram.get(str(v), 0)+1
    hits = sum(v >= threshold for v in values)
    return dict(reads=n, hits=hits, fraction=hits/n if n else None,
                fraction_at_least={str(t): sum(v >= t for v in values)/n for t in grid} if n else {},
                quantiles={f'p{q}': nearest(q) for q in (0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100)} if n else {},
                mean=sum(values)/n if n else None, longest_match_histogram=histogram)


def rank(rows, min_bases):
    """Highest fraction first; equal fractions break by usable bases, then accession."""
    exact = {r['run']: Fraction(r['screen']['hits'], r['screen']['reads']) if r['screen']['reads'] else Fraction(-1) for r in rows}
    ordered = sorted(rows, key=lambda r: (-exact[r['run']], -r['usable_bases'], r['run']))
    ranking = []
    for i, r in enumerate(ordered, 1):
        ranking.append(dict(rank=i, run=r['run'], fraction=r['screen']['fraction'], hits=r['screen']['hits'],
                            reads=r['screen']['reads'], usable_bases=r['usable_bases'],
                            eligible=r['usable_bases'] >= min_bases and r['screen']['reads'] > 0,
                            tied_with=[o['run'] for o in ordered if o is not r and exact[o['run']] == exact[r['run']]]))
    eligible = [r for r in ranking if r['eligible']]
    if not eligible:raise ValueError('No run has enough usable bases')
    selected = eligible[0]
    ties = [r['run'] for r in eligible[1:] if exact[r['run']] == exact[selected['run']]]
    return ranking, selected, ties


def save_new(path, value):
    path = Path(path);tmp = path.with_name(path.name+'.tmp')
    with tmp.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False);handle.write('\n');handle.flush();os.fsync(handle.fileno())
    try:os.link(tmp, path)
    finally:tmp.unlink()


def screen(manifest, output, work, *, root, ropebwt3, index, threads=16, pairs=PAIRS,
           threshold=THRESHOLD, min_bases=MIN_BASES, seed=SEED, expected_runs=77, runner=default_runner):
    output, work = Path(output), Path(work)
    if output.exists():raise FileExistsError(output)
    if isinstance(threshold, bool) or not isinstance(threshold, int) or threshold < 1:raise ValueError('Integer threshold required')
    if not 1 <= threads <= 48:raise ValueError('Bounded thread count required')
    manifest_hash = sha256(manifest)
    deposited = runs(manifest, root)
    if expected_runs is not None and len(deposited) != expected_runs:raise ValueError('Unexpected skin run count')
    for (path, expected), label in ((ropebwt3, 'ropebwt3'), (index, 'index')):
        if sha256(path) != expected:raise ValueError('Pinned checksum differs: '+label)
    work.mkdir(parents=True, exist_ok=False)
    started = time.time()
    jobs = [(run, r1, r2, pairs, seed) for run, (r1, r2) in deposited.items()]
    with multiprocessing.get_context('spawn').Pool(min(threads, len(jobs))) as pool:
        sampled = pool.map(sample_run, jobs)
    owners = []
    with (work/'screen.fa').open('w') as out:
        for r, (row, sample) in enumerate(sampled):
            for p, (_, _, a, b) in enumerate(sample):
                for mate, sequence in enumerate((a, b)):
                    for _, _, part in segments(sequence):
                        out.write(f'>q{len(owners)}\n{part}\n');owners.append((r, 2*p+mate))
    if not owners:raise ValueError('Screen sample is empty')
    argv = [str(ropebwt3[0]), 'mem', '-t', str(threads), '-l1', str(index[0]), str(work/'screen.fa')]
    match_started = time.time()
    runner(argv, work/'hits.tsv', work/'ropebwt3.stderr')
    match_seconds = time.time()-match_started
    best = longest_matches(work/'hits.tsv', len(owners))
    per_read = [[0]*(2*row['sampled_pairs']) for row, _ in sampled]
    for (r, i), value in zip(owners, best):
        per_read[r][i] = max(per_read[r][i], value)
    rows = [dict(row, sampled_reads=2*row['sampled_pairs'], screen=summarize(values, threshold))
            for (row, _), values in zip(sampled, per_read)]
    ranking, selected, ties = rank(rows, min_bases)
    result = dict(status='screened', cohort='guccione_skin', reference='human472', publishable_timing=False, note=CAVEAT,
                  parameters=dict(seed=seed, pairs_per_run=pairs, reads_per_run=2*pairs, selection='bottom-k pairs by SHA256(seed:run:pair-name), both mates as independent reads; every pair when a run has fewer',
                                  threshold=threshold, threshold_justification=JUSTIFICATION, grid=list(GRID), min_usable_bases=min_bases,
                                  usable_bases='ACGT bases after uppercasing; other bases split reads as in the historical cohorts',
                                  ranking='fraction of sampled reads with longest exact match >= threshold, descending; ties by usable bases descending, then run accession'),
                  manifest=str(manifest), manifest_sha256=manifest_hash, reads_root=str(root),
                  ropebwt3=dict(path=str(ropebwt3[0]), sha256=ropebwt3[1]), index=dict(path=str(index[0]), sha256=index[1]),
                  command=argv, threads=threads, host=socket.gethostname(), job_id=os.environ.get('SLURM_JOB_ID'),
                  work=str(work), screen_fasta_sha256=sha256(work/'screen.fa'), segments=len(owners),
                  started=started, match_seconds_untimed=match_seconds, finished=time.time(),
                  runs=rows, ranking=ranking, selected=selected['run'], selected_ties=ties,
                  tie_note='Selected run shares its exact fraction with: '+', '.join(ties) if ties else 'No eligible run shares the selected fraction')
    (work/'hits.tsv').unlink()
    save_new(output, result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--root', type=Path, required=True, help='the download root of data/download/download_skin.py')
    parser.add_argument('--manifest', type=Path, help='ROOT/skin-verified.json by default')
    parser.add_argument('--ropebwt3', type=Path, required=True)
    parser.add_argument('--fmd', type=Path, required=True, help='human472.fmd')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--pairs', type=int, default=PAIRS)
    parser.add_argument('--threshold', type=int, default=THRESHOLD)
    x = parser.parse_args(argv)
    r = screen(x.manifest or x.root/'skin-verified.json', x.output, x.work, root=x.root,
               ropebwt3=(x.ropebwt3, ROPEBWT3_SHA256), index=(x.fmd, INDEX_SHA256), threads=x.threads,
               pairs=x.pairs, threshold=x.threshold)
    print(json.dumps(dict(selected=r['selected'], ties=r['selected_ties'], top=r['ranking'][:5]), indent=2), flush=True)


if __name__ == '__main__':
    main()
