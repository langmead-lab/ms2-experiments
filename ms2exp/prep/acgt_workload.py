#!/usr/bin/env python3
"""ACGT-only whole-read workloads cut from a read pool, with sliced truth (the explore-acgt-60mb-v1 sets).

The PacBio and Nanopore 60 Mbase sets of the freeze are d13_prep.py caps of
campaign-20261012/inputs/COHORT/explore-acgt-60mb-v1, which this module made from the pools of
ms2exp.prep.pool and their truth (ms2exp.prep.truth_pool). Their manifests are in
data/evidence/rockfish/results/campaign-20261012/inputs/. Like the pools, this is provenance and is not
rerun by data/prep/run_d13.sh.

A workload is the first records of a verified source, in source order, until at least TARGET bases are
taken. A record qualifies when it is a whole original read made only of A, C, G and T after uppercasing;
nothing is split or trimmed and there is no minimum length. Records are renumbered q0, q1, ...;
source.tsv maps them back. Truth is sliced from the source's truth, since a read's matching statistics do
not depend on other reads.

Sources:
  core COHORT: POOLS/COHORT/core.fa, the seeded SeqKit sample of the cohort's run, and its truth
    TRUTH/COHORT/core.u32. A record is whole when its core.map.tsv row starts at 0 and ends at the read's
    length in original-read-names.tsv. When CELLS/core/COHORT/TOOL/1/result.json exists for the ms2 and
    TeraMS tools of the real-hprc2 campaign, their whole-output checks against this truth are recorded.
  workload DIR: a workload already filtered to whole ACGT-only reads.

With --all-if-short, a source holding fewer qualifying bases than the target gives all of them, and the
manifest records the shortfall; without it, that is an error.

Usage:
  python -m ms2exp.prep.acgt_workload core COHORT TARGET_BASES DESTINATION --root REAL [--all-if-short]
  python -m ms2exp.prep.acgt_workload workload SOURCE_DIR TARGET_BASES DESTINATION [--all-if-short]
REAL holds inputs/ (the pools), truth/ and optionally cells/, as real-hprc2-1d580b6 on Rockfish does.
"""
import array
import json
from pathlib import Path
import sys
from ms2exp.run.common import atomic_save
from ms2exp.sim.chr22_sim import sha256

FILTER = 'acgt-only-whole-reads'
# The real-hprc2 campaign directory on Rockfish; --root replaces it.
REAL = Path('/scratch16/blangme2/ms-bakery/results/real-hprc2-1d580b6')
CORE_CELLS = ('Orbit-ms', 'Orbit-tms-psi', 'Orbit-tms-dual', 'TeraMS-psi', 'TeraMS-dual')


def fasta(path):
    """(name, sequence) records, sequences joined across lines."""
    name, chunks = None, []
    with Path(path).open() as stream:
        for line in stream:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, ''.join(chunks)


def read_offsets(path):
    bounds = array.array('Q', Path(path).read_bytes())
    if sys.byteorder != 'little':
        bounds.byteswap()
    return bounds


def select(records, bounds, whole, target, all_if_short=False):
    """Indices and sequences of qualifying source records until at least target bases.

records yields (name, sequence) in source order with names q0, q1, ...; whole(i) says
whether source record i is a whole original read. With all_if_short, a source with
fewer qualifying bases gives all of them instead of an error.
"""
    kept, skipped, bases = [], dict(records=0, bases=0, not_whole=0, non_acgt=0), 0
    for i, (name, sequence) in enumerate(records):
        if bases >= target:
            break
        if name != f'q{i}' or len(sequence) != bounds[i + 1] - bounds[i]:
            raise ValueError(f'Source record {i} differs from its boundaries')
        upper = sequence.upper()
        reason = 'not_whole' if not whole(i) else 'non_acgt' if upper.strip('ACGT') else None
        if reason:
            skipped[reason] += 1
            skipped['records'] += 1
            skipped['bases'] += len(sequence)
            continue
        kept.append((i, upper))
        bases += len(upper)
    if bases < target and not all_if_short:
        raise ValueError('Source holds fewer qualifying bases than the target')
    return kept, skipped


def write(destination, kept, bounds, truth):
    """reads.fa, one.fa, offsets.u64, truth.u32 (sliced) and source.tsv; returns file records."""
    destination.mkdir(parents=True, exist_ok=False)
    prefix = array.array('Q', [0])
    with (destination/'reads.fa').open('x') as out:
        for j, (_, sequence) in enumerate(kept):
            out.write(f'>q{j}\n{sequence}\n')
            prefix.append(prefix[-1] + len(sequence))
    (destination/'one.fa').write_text(f'>q0\n{kept[0][1]}\n')
    if sys.byteorder != 'little':
        prefix.byteswap()
    (destination/'offsets.u64').write_bytes(prefix.tobytes())
    with Path(truth).open('rb') as source, (destination/'truth.u32').open('xb') as out:
        for i, sequence in kept:
            source.seek(4 * bounds[i])
            block = source.read(4 * len(sequence))
            if len(block) != 4 * len(sequence):
                raise ValueError('Source truth truncated')
            values = array.array('I', block)
            if sys.byteorder != 'little':
                values.byteswap()
            if any(x > len(sequence) - k for k, x in enumerate(values)):
                raise ValueError('Truth value exceeds the rest of its read')
            out.write(block)
    (destination/'source.tsv').write_text(''.join(f'q{j}\tq{i}\n' for j, (i, _) in enumerate(kept)))
    names = ('reads.fa', 'one.fa', 'truth.u32', 'offsets.u64', 'source.tsv')
    return {name: dict(bytes=(destination/name).stat().st_size, sha256=sha256(destination/name)) for name in names}


def core(cohort, target, destination, all_if_short=False, root=None):
    root = Path(root or REAL)
    inputs, truth_dir = root/'inputs'/cohort, root/'truth'/cohort
    complete = json.loads((inputs/'complete.json').read_text())
    meta = json.loads((truth_dir/'core.json').read_text())
    expected = complete['sets']['core']
    if not meta.get('passed') or meta['input']['sha256'] != expected['sha256']:
        raise ValueError('Core truth did not pass or describes other reads')
    if sha256(inputs/'core.fa') != expected['sha256'] or sha256(inputs/'core.map.tsv') != expected['map_sha256']:
        raise ValueError('Core reads or map differ from their record')
    bounds = read_offsets(inputs/'core.offsets.u64')
    if len(bounds) != expected['reads'] + 1 or bounds[-1] != expected['bases']:
        raise ValueError('Core boundaries differ')
    if (truth_dir/'core.u32').stat().st_size != 4 * expected['bases']:
        raise ValueError('Core truth size differs')
    lengths = {}
    for line in (inputs/'original-read-names.tsv').read_text().splitlines():
        fields = line.split('\t')
        lengths[fields[0]] = int(fields[-1])
    rows = [line.split('\t')[:4] for line in (inputs/'core.map.tsv').read_text().splitlines()]
    keys = [row[1] for row in rows]
    if len(set(keys)) != len(keys):
        raise ValueError('An original read has several core records')
    def whole(i):
        q, key, start, end = rows[i]
        return q == f'q{i}' and int(start) == 0 and int(end) == lengths[key]
    kept, skipped = select(fasta(inputs/'core.fa'), bounds, whole, target, all_if_short)
    files = write(destination, kept, bounds, truth_dir/'core.u32')
    cells = {}
    for name in CORE_CELLS:
        path = root/'cells/core'/cohort/name/'1/result.json'
        if path.exists():
            result = json.loads(path.read_text())
            checks = [p['raw']['full']['check'] for p in result.get('pairs', [])]
            cells[name] = dict(path=str(path), status=result.get('status'),
                               full_checks_failed=[c.get('failed') for c in checks],
                               positions=[c.get('positions') for c in checks])
    truth = dict(source=str(truth_dir/'core.u32'), source_sha256=sha256(truth_dir/'core.u32'),
                 metadata=str(truth_dir/'core.json'), metadata_sha256=sha256(truth_dir/'core.json'),
                 generator=meta['command']['argv'], generator_check=meta['command']['check'],
                 method='ropebwt3 mem -l1 SMEMs on the Human472 FMD converted to MS; sliced per read',
                 cross_tool_validation=cells)
    source = dict(kind='core', path=str(inputs/'core.fa'), sha256=expected['sha256'],
                  map=str(inputs/'core.map.tsv'), map_sha256=expected['map_sha256'],
                  sampling=complete['commands'][:2], reads=expected['reads'], bases=expected['bases'])
    return finish(destination, cohort, target, kept, skipped, files, source, truth, run=None)


def workload(source_dir, target, destination, all_if_short=False):
    source_dir = Path(source_dir)
    manifest = json.loads((source_dir/'manifest.json').read_text())
    if manifest.get('read_filter') != FILTER or manifest.get('reference') != 'human472':
        raise ValueError('An ACGT-only Human472 workload is required')
    for name in ('reads.fa', 'truth.u32', 'offsets.u64'):
        entry = manifest['files'][name]
        if (source_dir/name).stat().st_size != entry['bytes'] or sha256(source_dir/name) != entry['sha256']:
            raise ValueError('Source checksum or size differs: ' + name)
    bounds = read_offsets(source_dir/'offsets.u64')
    kept, skipped = select(fasta(source_dir/'reads.fa'), bounds, lambda i: True, target, all_if_short)
    files = write(destination, kept, bounds, source_dir/'truth.u32')
    source = dict(kind='workload', path=str(source_dir), manifest_sha256=sha256(source_dir/'manifest.json'),
                  reads=manifest['reads'], bases=manifest['bases'], layout=manifest.get('layout'))
    truth = dict(source=str(source_dir/'truth.u32'), source_sha256=manifest['files']['truth.u32']['sha256'],
                 method='sliced per read from the source workload truth (' + manifest.get('truth', 'recorded there') + ')')
    return finish(destination, manifest['cohort'], target, kept, skipped, files, source, truth, run=manifest.get('run'))


def finish(destination, cohort, target, kept, skipped, files, source, truth, run):
    if not kept:
        raise ValueError('No qualifying records')
    lengths = sorted(len(s) for _, s in kept)
    result = dict(status='prepared', reference='human472', cohort=cohort, mode='exploratory',
                  reads=len(kept), bases=sum(lengths), one_bases=len(kept[0][1]), target_bases=target,
                  read_lengths=dict(min=lengths[0], median=lengths[len(lengths) // 2], max=lengths[-1]),
                  read_filter=FILTER, min_record_length=None,
                  selection=('first records of the source in its fixed order whose original read is whole and '
                             'consists only of A, C, G and T after uppercasing, until the target is reached; '
                             'nothing split or trimmed; no minimum length; renumbered q0, q1, ...'),
                  shortfall=(dict(target_bases=target, bases=sum(lengths), rule='all qualifying records taken')
                             if sum(lengths) < target else None),
                  skipped=skipped, source=source, truth=truth, run=run, reuse_old_timings=False, files=files)
    result = {k: v for k, v in result.items() if v is not None or k == 'min_record_length'}
    atomic_save(destination/'manifest.json', result)
    return result


def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('kind', choices=('core', 'workload'))
    p.add_argument('source', help='COHORT for core, SOURCE_DIR for workload')
    p.add_argument('target', type=int)
    p.add_argument('destination', type=Path)
    p.add_argument('--root', help='the pool campaign directory (core only)')
    p.add_argument('--all-if-short', action='store_true')
    a = p.parse_args(argv)
    if a.kind == 'core':
        out = core(a.source, a.target, a.destination, a.all_if_short, root=a.root)
    else:
        out = workload(a.source, a.target, a.destination, a.all_if_short)
    print(json.dumps({k: v for k, v in out.items() if k not in ('files', 'truth', 'source')}, indent=2))


if __name__ == '__main__':
    main()
