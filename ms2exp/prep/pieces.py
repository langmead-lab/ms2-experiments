#!/usr/bin/env python3
"""Genome query workloads and their 100 kbase piece sets (Fig ops genome queries), with truth.

The genome query sets of Fig ops were made in three steps, then sampled by d13_prep.py sample:
  workload: split every record of a held-out genome at runs of non-ACGT characters (uppercased first)
    into q-named records (reads.fa, offsets.u64, map.tsv, reads.json);
  truth: ropebwt3 mem -l1 on the index's FMD over the whole workload, converted to per-position
    matching statistics by ms2exp.validate.ms.check and bounds-checked (truth.u32, truth.json,
    manifest.json);
  pieces: cut every record into consecutive pieces of at most --length bases (the last piece of a
    record may be shorter), with truth min(MS[i], end - i), where end is the piece's end in its record.
    Every prefix of a string that occurs in the index also occurs, so this is the pieces' exact truth.
The Fig ops PacBio sets used workload (with the simulated reads as the source, cohort fig_ops2) and
truth. The frozen records are in data/evidence/rockfish/results/campaign-20261012/inputs/fig_ops*/.
The pieces step needs numpy.

Usage:
  python -m ms2exp.prep.pieces workload --reference R --name NAME --source FASTA --destination DIR
                                        [--cohort fig_ops] [--source-note TEXT] [--max-record-length N]
  python -m ms2exp.prep.pieces truth --workload DIR --ropebwt3 PATH --fmd PATH --work DIR [--threads N]
  python -m ms2exp.prep.pieces pieces --workload DIR --length 100000 --destination DIR
"""
import argparse
import array
import json
import os
from pathlib import Path
import re
import shutil
import sys
import time

from ms2exp.run import common
from ms2exp.sim.chr22_sim import sha256
from ms2exp.validate.ms import check

FILTER = 'acgt-split-at-non-acgt'
WORKLOAD_FILES = ('reads.fa', 'truth.u32', 'offsets.u64')
# The held-out sources of the two references and their notes as recorded in the frozen manifests.
REFERENCES = {
    'mtb152': dict(source_note=('TeraTools paper MS benchmark, Nate Brown 95/5 split: the 5 percent held out as '
                                'patterns (8 genomes, each followed by its reverse complement)')),
    'human472': dict(source_note=('Arab Pangenome Reference APR008 haplotype 1, GenBank GCA_050494745.1 '
                                  '(apr008.1 v1), via inputs/heldout/apr-APR008-hap1 (manifest.json there has '
                                  'provenance and held-out evidence)')),
}


def save_new(path, value):
    """Write JSON to a new file, refusing to replace one."""
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp')
    with tmp.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    try:
        os.link(tmp, path)
    finally:
        tmp.unlink()


def read_u64(path):
    values = array.array('Q')
    values.frombytes(Path(path).read_bytes())
    if sys.byteorder != 'little':
        values.byteswap()
    return values


def write_u64(path, values):
    a = array.array('Q', values)
    if sys.byteorder != 'little':
        a.byteswap()
    Path(path).write_bytes(a.tobytes())


# Workloads, on the dev node: ACGT-only records named q0, q1, ...

def fasta_pieces(path, block=1 << 24):
    """Yield ('header', name) and ('seq', bytes) items, reading lines in bounded blocks."""
    at_line_start = True
    with Path(path).open('rb') as f:
        while True:
            part = f.readline(block)
            if not part:
                return
            ends = part.endswith(b'\n')
            if at_line_start and part.startswith(b'>'):
                header = part
                while not header.endswith(b'\n'):
                    more = f.readline(block)
                    if not more:
                        break
                    header += more
                fields = header[1:].split()
                yield 'header', (fields[0].decode() if fields else '')
                at_line_start = True
                continue
            yield 'seq', part.rstrip(b'\r\n')
            at_line_start = ends


def make_workload(reference, name, source, destination, max_record_length=None, block=1 << 24,
                  cohort='fig_ops', source_note=None):
    """Split every record at runs of non-ACGT characters (uppercased first) into q-named records.

With max_record_length, pieces longer than that are cut into consecutive records of at
most that length. Writes reads.fa (two-line FASTA), offsets.u64, map.tsv (q name,
source record, 0-based start and exclusive end in it) and reads.json. Fig Ops2 workloads
(cohort fig_ops2) take simulated read sets as the source, which are ACGT-only already, so
each read becomes one record of the same name; source_note then describes the read set.
"""
    source = Path(source).resolve(strict=True)
    destination = Path(destination)
    if max_record_length is not None and max_record_length < 1:
        raise ValueError('max_record_length must be positive')
    destination.mkdir(parents=True, exist_ok=False)
    bad = re.compile(rb'[^ACGT]+')
    offsets, total, dropped, records_in = [0], 0, 0, 0
    with (destination/'reads.fa').open('xb') as out, (destination/'map.tsv').open('x') as mapping:
        state = dict(record=None, pos=0, open=False, start=0, length=0)

        def close():
            if state['open']:
                out.write(b'\n')
                q = len(offsets) - 1
                mapping.write(f"q{q}\t{state['record']}\t{state['start']}\t{state['start'] + state['length']}\n")
                offsets.append(offsets[-1] + state['length'])
                state['open'] = False

        def emit(data, start):
            """Append ACGT bytes that begin at source position start, opening and cutting records as needed."""
            while data:
                if not state['open']:
                    out.write(f'>q{len(offsets) - 1}\n'.encode())
                    state.update(open=True, start=start, length=0)
                room = len(data) if max_record_length is None else max_record_length - state['length']
                take = data[:room]
                out.write(take)
                state['length'] += len(take)
                start += len(take)
                data = data[len(take):]
                if max_record_length is not None and state['length'] == max_record_length:
                    close()

        for kind, value in fasta_pieces(source, block):
            if kind == 'header':
                close()
                records_in += 1
                state.update(record=value, pos=0)
                continue
            if state['record'] is None:
                raise ValueError('Sequence before the first header')
            chunk = value.upper()
            last = 0
            for m in bad.finditer(chunk):
                emit(chunk[last:m.start()], state['pos'] + last)
                close()
                dropped += m.end() - m.start()
                last = m.end()
            emit(chunk[last:], state['pos'] + last)
            state['pos'] += len(chunk)
            total += len(chunk)
        close()
    write_u64(destination/'offsets.u64', offsets)
    files = {n: dict(bytes=(destination/n).stat().st_size, sha256=sha256(destination/n))
             for n in ('reads.fa', 'offsets.u64', 'map.tsv')}
    lengths = [b - a for a, b in zip(offsets, offsets[1:])]
    record = dict(status='reads-prepared', reference=reference, cohort=cohort, set=name, read_filter=FILTER,
                  source=str(source), source_sha256=sha256(source), source_note=source_note or REFERENCES.get(reference, {}).get('source_note'),
                  rule=('uppercase; split every record at each run of characters other than A, C, G and T, which are '
                        'dropped; empty pieces dropped' + ('' if max_record_length is None else
                                                          f'; pieces cut into records of at most {max_record_length} bases')
                        + '; records renamed q0, q1, ... in source order (map.tsv gives each source record and '
                          '0-based half-open range)'),
                  max_record_length=max_record_length, source_records=records_in, source_bases=total,
                  dropped_non_acgt=dropped, reads=len(offsets) - 1, bases=offsets[-1],
                  longest=max(lengths) if lengths else 0, shortest=min(lengths) if lengths else 0,
                  created=time.time(), files=files)
    if record['bases'] + dropped != total:
        raise ValueError('Base accounting differs')
    save_new(destination/'reads.json', record)
    return record


# Truth, in a Slurm job: ropebwt3 -l1 on the reference's FMD, converted by ms2exp.validate.ms.check.

def bounded(truth, offsets):
    """Every ACGT base occurs in the reference and no match runs past its read (streamed in blocks)."""
    block = 1 << 24
    with Path(truth).open('rb') as f:
        for i in range(len(offsets) - 1):
            n = offsets[i + 1] - offsets[i]
            j = 0
            while j < n:
                count = min(block, n - j)
                ms = array.array('I')
                ms.frombytes(f.read(4 * count))
                if sys.byteorder != 'little':
                    ms.byteswap()
                if len(ms) != count or any(not 1 <= v <= n - (j + k) for k, v in enumerate(ms)):
                    raise ValueError(f'MS truth out of range for q{i}')
                j += count


def truth(workload, rope, fmd, work, threads=16, batch='10g', timeout=12 * 3600, check=check):
    """Write truth.u32, truth.json and manifest.json into a prepared workload directory.

rope and fmd are (path, sha256) pairs; both checksums are checked first."""
    workload = Path(workload).resolve(strict=True)
    reads = json.loads((workload/'reads.json').read_text())
    if reads.get('status') != 'reads-prepared' or reads.get('read_filter') != FILTER:
        raise ValueError('Prepared Fig Ops reads required')
    for name in ('truth.u32', 'truth.json', 'manifest.json'):
        if (workload/name).exists():
            raise FileExistsError(workload/name)
    for name, entry in reads['files'].items():
        if sha256(workload/name) != entry['sha256']:
            raise ValueError('Prepared file differs: ' + name)
    for (path, expected), label in ((rope, 'ropebwt3'), (fmd, 'FMD')):
        if sha256(path) != expected:
            raise ValueError('Pinned checksum differs: ' + label)
    offsets = read_u64(workload/'offsets.u64')
    if len(offsets) != reads['reads'] + 1 or offsets[-1] != reads['bases']:
        raise ValueError('Read boundaries differ from reads.json')
    work = Path(work)/f"figops-{workload.name}-{os.environ.get('SLURM_JOB_ID', 'local')}"
    work.mkdir(parents=True, exist_ok=False)
    started = time.time()
    # Every SMEM of length 1 or more; one batch (-K) holds the whole genome so threads share it.
    argv = [rope[0], 'mem', '-t', threads, '-l1', '-K', batch, fmd[0], workload/'reads.fa']
    row = common.execute(argv, work, 'truth', timeout)
    partial = work/'truth.u32'
    row['check'] = check('ropebwt3', work, row, partial, offsets, write=True)
    bounded(partial, offsets)
    row['stdout_sha256'] = sha256(row['stdout'])
    row['stdout_bytes'] = Path(row['stdout']).stat().st_size
    Path(row['stdout']).unlink()
    row['raw_output_retained'] = False
    common.save(work/'truth.command.json', row)
    with partial.open('rb') as src, (workload/'truth.u32').open('xb') as dst:
        shutil.copyfileobj(src, dst, 16 << 20)
        dst.flush()
        os.fsync(dst.fileno())
    digest = sha256(workload/'truth.u32')
    if digest != sha256(partial):
        raise ValueError('Truth copy differs')
    partial.unlink()
    verdict = dict(passed=True, command=row, reference=reads['reference'],
                   reads_manifest_sha256=sha256(workload/'reads.json'), sha256=digest,
                   ropebwt3=dict(path=str(rope[0]), sha256=rope[1]), fmd=dict(path=str(fmd[0]), sha256=fmd[1]),
                   derivation=('ropebwt3 mem -l1 SMEMs converted to per-position MS by ms2exp.validate.ms.check; '
                               'every value checked to satisfy 1 <= ms[j] <= length - j'),
                   job_id=os.environ.get('SLURM_JOB_ID'), started=started, finished=time.time())
    save_new(workload/'truth.json', verdict)
    names = ('reads.fa', 'offsets.u64', 'map.tsv', 'truth.u32', 'truth.json')
    manifest = dict(status='prepared', reference=reads['reference'], cohort=reads['cohort'], set=reads['set'],
                    reads=reads['reads'], bases=reads['bases'], read_filter=FILTER,
                    selection=f"{reads['source_note']}; {reads['rule']}", source=reads['source'],
                    reads_manifest_sha256=sha256(workload/'reads.json'),
                    truth=(f"ropebwt3 -l1 (binary {rope[1][:12]}) on the {reads['reference']} FMD "
                           f"({fmd[1][:12]}), job {os.environ.get('SLURM_JOB_ID')}"),
                    truth_metadata_sha256=sha256(workload/'truth.json'), reuse_old_timings=False,
                    files={n: dict(bytes=(workload/n).stat().st_size, sha256=sha256(workload/n)) for n in names})
    save_new(workload/'manifest.json', manifest)
    return manifest


def pieces(workload, length, destination, block=1 << 24):
    """Cut every record of a prepared workload into consecutive pieces of at most length bases.

Truth comes from the source workload's: a piece's MS at position i is min(MS[i], end - i),
where end is the piece's end in its record, because every prefix of a string that occurs
in the index also occurs. map.tsv gives each piece's source q name and 0-based half-open
range in it.
"""
    import numpy as np
    workload = Path(workload).resolve(strict=True)
    manifest = json.loads((workload/'manifest.json').read_text())
    if manifest.get('status') != 'prepared' or manifest.get('read_filter') != FILTER:
        raise ValueError('A prepared Fig Ops workload is required')
    for name in WORKLOAD_FILES:
        if sha256(workload/name) != manifest['files'][name]['sha256']:
            raise ValueError('Workload file differs from its manifest: ' + name)
    if length < 1:
        raise ValueError('length must be positive')
    offsets = read_u64(workload/'offsets.u64')
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    out_offsets, clamped = [0], 0
    with (destination/'reads.fa').open('xb') as out, (destination/'map.tsv').open('x') as mapping, \
            (destination/'truth.u32').open('xb') as truth_out, (workload/'truth.u32').open('rb') as truth_in:
        q, seq = -1, []

        def flush():
            nonlocal clamped
            data = b''.join(seq)
            if len(data) != offsets[q + 1] - offsets[q]:
                raise ValueError(f'Record q{q} length differs from offsets.u64')
            for start in range(0, len(data), length):
                piece = data[start:start + length]
                n = len(piece)
                ms = np.frombuffer(truth_in.read(4 * n), dtype='<u4')
                if len(ms) != n:
                    raise ValueError('Truth ends early')
                cut = np.minimum(ms, np.arange(n, 0, -1, dtype=np.uint32))
                if cut.min() < 1:
                    raise ValueError(f'MS truth below 1 in q{q}')
                clamped += int(np.count_nonzero(cut != ms))
                p = len(out_offsets) - 1
                out.write(f'>q{p}\n'.encode() + piece + b'\n')
                truth_out.write(cut.astype('<u4').tobytes())
                mapping.write(f'q{p}\tq{q}\t{start}\t{start + n}\n')
                out_offsets.append(out_offsets[-1] + n)

        for kind, value in fasta_pieces(workload/'reads.fa', block):
            if kind == 'header':
                if q >= 0:
                    flush()
                q, seq = q + 1, []
                if value != f'q{q}':
                    raise ValueError('Workload records are not q0, q1, ... in order')
            else:
                seq.append(value)
        if q >= 0:
            flush()
        if q + 1 != len(offsets) - 1 or truth_in.read(1):
            raise ValueError('Record count or truth length differs from the workload')
    write_u64(destination/'offsets.u64', out_offsets)
    bounded(destination/'truth.u32', out_offsets)
    names = ('reads.fa', 'offsets.u64', 'map.tsv', 'truth.u32')
    lengths = [b - a for a, b in zip(out_offsets, out_offsets[1:])]
    derivation = (f"truth of {workload.name} (manifest {sha256(workload/'manifest.json')[:12]}) clamped per piece: "
                  'min(MS[i], end - i)')
    save_new(destination/'truth.json', dict(passed=True, derivation=derivation, source_workload=str(workload),
                                            source_truth=manifest.get('truth'), clamped_positions=clamped,
                                            created=time.time()))
    record = dict(status='prepared', reference=manifest['reference'], cohort=manifest.get('cohort', 'fig_ops'),
                  set=destination.name, reads=len(out_offsets) - 1, bases=out_offsets[-1], read_filter=FILTER,
                  selection=(f"{manifest.get('selection')}; then every record cut into consecutive pieces of at most "
                             f'{length} bases, renamed q0, q1, ... in order (map.tsv gives each source record and '
                             '0-based half-open range)'),
                  source=str(workload), piece_length=length, longest=max(lengths), shortest=min(lengths),
                  truth=derivation, truth_metadata_sha256=sha256(destination/'truth.json'), reuse_old_timings=False,
                  files={n: dict(bytes=(destination/n).stat().st_size, sha256=sha256(destination/n)) for n in names})
    if record['bases'] != manifest['bases']:
        raise ValueError('Base count differs from the source workload')
    save_new(destination/'manifest.json', record)
    return record


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='action', required=True)
    w = sub.add_parser('workload')
    w.add_argument('--reference', required=True)
    w.add_argument('--name', required=True)
    w.add_argument('--source', required=True)
    w.add_argument('--destination', required=True)
    w.add_argument('--cohort', default='fig_ops')
    w.add_argument('--source-note')
    w.add_argument('--max-record-length', type=int)
    t = sub.add_parser('truth')
    t.add_argument('--workload', required=True)
    t.add_argument('--ropebwt3', required=True)
    t.add_argument('--fmd', required=True)
    t.add_argument('--work', required=True, help='directory for the ropebwt3 run')
    t.add_argument('--threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 16)))
    c = sub.add_parser('pieces')
    c.add_argument('--workload', required=True)
    c.add_argument('--length', type=int, required=True)
    c.add_argument('--destination', required=True)
    a = p.parse_args(argv)
    if a.action == 'workload':
        m = make_workload(a.reference, a.name, a.source, a.destination, a.max_record_length, cohort=a.cohort,
                          source_note=a.source_note)
    elif a.action == 'truth':
        m = truth(a.workload, (a.ropebwt3, sha256(a.ropebwt3)), (a.fmd, sha256(a.fmd)), a.work, threads=a.threads)
    else:
        m = pieces(a.workload, a.length, a.destination)
    print(json.dumps({k: v for k, v in m.items() if k != 'files'}), flush=True)


if __name__ == '__main__':
    main()
