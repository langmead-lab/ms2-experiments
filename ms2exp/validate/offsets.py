"""Locate (SMEM offset) validation and canonical SMEM location files.

SMEM means super-maximal exact match. A canonical file (format canonical-smem-v1) is UTF-8 text
compressed with gzip (mtime 0, no stored file name). It has one line per SMEM, tab-separated as
ropebwt3 mem -p N writes it:

    read  start  end  count  n_pos  contig:strand:offset ...

Reads appear in workload order (q0, q1, ...), a read's lines by start, and the positions of a line
sorted by contig name (as a string), strand and offset. start and end are 0-based and half-open in
the read; offset is the 0-based forward-strand start of the match, as ropebwt3 reports it for
either strand. n_pos is min(cap, count). A line with count 0 has no n_pos or positions and marks a
read character that occurs nowhere in the reference. The SHA256 of the uncompressed text
identifies the content; the compressed bytes are only deterministic for one zlib build, so they
are never compared across tools.

check_smems() applies the cross-tool policy (cross-tool-v1) to output in ropebwt3 mem -p N form,
which both ropebwt3 and ms2 query --report smem-all write. It needs no reference sequence: it
checks structure, the SMEM intervals and the matching statistics they imply against MS truth, and
that every position lies inside its contig, and writes the canonical form. Counts and positions
are certified by exact agreement between the two tools (compare_files, run by ms2exp.collect);
positions of SMEMs whose count exceeds the cap may differ between tools and are not compared.
"""
import array
import gzip
import hashlib
import os
from pathlib import Path
import sys

FORMAT = 'canonical-smem-v1'


def canonical_line(name, start, end, count, positions=()):
    fields = [name, str(start), str(end), str(count)]
    if count:
        fields.append(str(len(positions)))
        fields.extend(f'{c}:{s}:{o}' for c, s, o in sorted(positions))
    return '\t'.join(fields) + '\n'


class CanonicalSink:
    """Hash canonical lines and, with a path, write them gzip-compressed.

The caller supplies lines already in canonical order.
"""
    def __init__(self, cap, path=None):
        self.cap, self.path = cap, Path(path) if path else None
        self.content, self.over = hashlib.sha256(), hashlib.sha256()
        self.lines = self.bytes = self.over_lines = self.over_locations = 0
        self.raw = self.stream = None
        if self.path:
            self.raw = self.path.open('xb')
            self.stream = gzip.GzipFile(filename='', mode='wb', fileobj=self.raw, mtime=0)

    def write(self, line, count, locations):
        data = line.encode()
        self.content.update(data)
        self.lines += 1
        self.bytes += len(data)
        if count > self.cap:
            self.over.update(data)
            self.over_lines += 1
            self.over_locations += locations
        if self.stream:
            self.stream.write(data)

    def close(self):
        summary = dict(format=FORMAT, sha256=self.content.hexdigest(), lines=self.lines,
                       bytes=self.bytes, cap=self.cap, over_cap_lines=self.over_lines,
                       over_cap_locations=self.over_locations, over_cap_sha256=self.over.hexdigest())
        if self.stream:
            self.stream.close()
            self.raw.flush()
            os.fsync(self.raw.fileno())
            self.raw.close()
            summary.update(compressed_bytes=self.path.stat().st_size, compressed_sha256=file_sha256(self.path))
        return summary

    def discard(self):
        if self.stream:
            self.stream.close()
            self.raw.close()
            self.path.unlink(missing_ok=True)


def file_sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(16 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read_canonical(path, sha256=None):
    """Yield (read, start, end, count, positions) and check the content SHA256 at the end."""
    h = hashlib.sha256()
    with gzip.open(path, 'rt', newline='\n') as stream:
        for line in stream:
            h.update(line.encode())
            fields = line.rstrip('\n').split('\t')
            start, end, count = map(int, fields[1:4])
            positions = ()
            if count:
                positions = tuple((c, s, int(o)) for c, s, o in (p.rsplit(':', 2) for p in fields[5:]))
                if int(fields[4]) != len(positions):
                    raise ValueError('Canonical position count differs')
            elif len(fields) != 4:
                raise ValueError('Malformed canonical zero-count line')
            yield fields[0], start, end, count, positions
    if sha256 is not None and h.hexdigest() != sha256:
        raise ValueError('Canonical output content differs from its recorded SHA256')


def read_rank(name):
    if not name.startswith('q') or not name[1:].isdigit():
        raise ValueError('Workload read identifier required')
    return int(name[1:])


def compare_files(orbit, rope, cap, examples=10):
    """Compare two canonical files given as (path, sha256); return counts and example differences.

Intervals and counts must agree for every SMEM, position sets for count <= cap.
Above the cap each tool may list different positions, which stay deferred.
"""
    left, right = read_canonical(*orbit), read_canonical(*rope)
    sentinel = (float('inf'),)
    def advance(stream):
        record = next(stream, None)
        return record, (sentinel if record is None else (read_rank(record[0]), record[1], record[2]))
    a, ka = advance(left)
    b, kb = advance(right)
    counts = dict(intervals=0, zero_count_lines=0, complete_sets=0, complete_set_locations=0,
                  over_cap_intervals=0, orbit_only=0, rope_only=0, count_differs=0,
                  listed_count_differs=0, position_set_differs=0)
    found = []
    def note(kind, record):
        counts[kind] += 1
        if len(found) < examples:
            found.append(dict(kind=kind, read=record[0], start=record[1], end=record[2]))
    while a is not None or b is not None:
        if ka < kb:
            note('orbit_only', a); a, ka = advance(left)
        elif kb < ka:
            note('rope_only', b); b, kb = advance(right)
        else:
            if a[3] != b[3]:
                note('count_differs', a)
            elif len(a[4]) != len(b[4]):
                note('listed_count_differs', a)
            elif a[3] == 0:
                counts['zero_count_lines'] += 1
            elif a[3] <= cap:
                if a[4] != b[4]:
                    note('position_set_differs', a)
                else:
                    counts['complete_sets'] += 1
                    counts['complete_set_locations'] += len(a[4])
            else:
                counts['over_cap_intervals'] += 1
            counts['intervals'] += a[3] > 0
            a, ka = advance(left)
            b, kb = advance(right)
    disagreements = sum(counts[k] for k in ('orbit_only', 'rope_only', 'count_differs',
                                            'listed_count_differs', 'position_set_differs'))
    return dict(status='agree' if disagreements == 0 else 'disagree', cap=cap,
                disagreements=disagreements, counts=counts, examples=found)


POLICY = 'cross-tool-v1'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_contigs(path):
    """Contig names and lengths, in index order, from a ropebwt3 .len.gz file."""
    contigs = {}
    with gzip.open(path, 'rt') as stream:
        for line in stream:
            fields = line.rstrip('\n').split('\t')
            require(len(fields) == 2 and fields[1].isdigit() and int(fields[1]) > 0
                    and fields[0] not in contigs, 'Malformed or duplicate sequence-length record')
            contigs[fields[0]] = int(fields[1])
    require(contigs, 'Empty sequence-length file')
    return contigs


def ms_from_smems(intervals, length):
    """Matching statistics implied by a read's SMEMs, given as sorted (start, end) pairs.

MS[i] is the largest end - i over SMEMs with start <= i < end, or 0: the longest
match starting at i, extended left as far as it occurs, is itself an SMEM.
"""
    ms, j, best = [0] * length, 0, None
    for i in range(length):
        while j < len(intervals) and intervals[j][0] <= i:
            best = intervals[j][1]
            j += 1
        if best is not None and best > i:
            ms[i] = best - i
    return ms


def check_smems(output, queries, ms_truth, contigs, cap, sink, tool):
    """Check ropebwt3-form SMEM output and pass its canonical lines to sink.

Reads must appear in workload order, each in one block of lines. A read with no
lines must have no SMEMs. Returns small evidence; counts and positions are
certified by cross-tool comparison of the canonical files.
"""
    # Production adapters accept caps 1 and 100; exploratory runs add 10.
    require(isinstance(cap, int) and not isinstance(cap, bool) and cap >= 1, 'Unsupported reporting cap')
    require(tool in ('orbit', 'rope'), 'Unsupported offset tool')
    names = list(queries)
    require(set(names) == set(ms_truth), 'Query/truth read IDs differ')
    rank = {name: i for i, name in enumerate(names)}
    totals = dict(intervals=0, zero=0, locations=0, ms=0, complete=0, complete_locations=0,
                  over=0, over_locations=0)
    state = dict(next=0)

    def close(name, hits):
        sequence, ms = queries[name], list(ms_truth[name])
        require(len(ms) == len(sequence), 'Truth length differs from query')
        require(all(0 <= x <= len(ms) - i for i, x in enumerate(ms)), 'Invalid MS truth')
        expected = {(i, i + x) for i, x in enumerate(ms) if x and (i == 0 or x >= ms[i - 1])}
        smems = sorted((s, e) for s, e, c, _ in hits if c)
        require(len(set((s, e) for s, e, _, _ in hits)) == len(hits), 'Duplicate interval')
        require(set(smems) == expected, 'Missing or extra SMEM intervals')
        require({s for s, e, c, _ in hits if not c} == {i for i, x in enumerate(ms) if not x},
                'Zero-count records differ from MS truth')
        require(ms_from_smems(smems, len(ms)) == ms, 'Matching statistics implied by SMEMs differ from truth')
        totals['ms'] += len(ms)
        for start, end, count, pos in sorted(hits):
            sink.write(canonical_line(name, start, end, count, pos), count, len(pos))

    def skip_to(index):
        # Reads between blocks wrote nothing, so they must have no SMEMs.
        while state['next'] < index:
            close(names[state['next']], [])
            state['next'] += 1

    current, hits = None, []
    with Path(output).open() as stream:
        for line in stream:
            fields = line.rstrip('\n').split('\t')
            name = fields[0]
            require(name in rank, 'Unexpected read ID')
            require(len(fields) >= 4, 'Malformed SMEM record')
            start, end, count = map(int, fields[1:4])
            if tool == 'rope' and count == 0 and start == 0 and end == 0 and len(fields) == 4:
                continue  # A zero-hit marker of ropebwt3 that records no interval.
            if name != current:
                if current is not None:
                    close(current, hits)
                    state['next'] = rank[current] + 1
                require(rank[name] >= state['next'], 'Reads out of workload order or split across blocks')
                skip_to(rank[name])
                current, hits = name, []
            require(0 <= start < end <= len(queries[name]), 'Invalid query interval')
            if count == 0:
                require(len(fields) == 4 and end == start + 1, 'Malformed zero-count record')
                hits.append((start, end, 0, ()))
                totals['zero'] += 1
                continue
            require(len(fields) >= 5 and count > 0, 'Malformed SMEM record')
            listed = int(fields[4])
            require(listed == min(cap, count) and len(fields) == 5 + listed, 'Incorrect number of capped locations')
            pos = []
            for item in fields[5:]:
                parts = item.rsplit(':', 2)
                require(len(parts) == 3 and parts[1] in ('+', '-') and parts[2].isdigit(),
                        'Position is not contig:strand:offset (an Orbit index without SEQ_IDS?)')
                contig, strand, offset = parts[0], parts[1], int(parts[2])
                require(contig in contigs and offset + end - start <= contigs[contig], 'Position outside its contig')
                pos.append((contig, strand, offset))
            require(len(set(pos)) == len(pos), 'Duplicate location')
            hits.append((start, end, count, tuple(pos)))
            totals['intervals'] += 1
            totals['locations'] += listed
            if count <= cap:
                totals['complete'] += 1
                totals['complete_locations'] += listed
            else:
                totals['over'] += 1
                totals['over_locations'] += listed
    if current is not None:
        close(current, hits)
        state['next'] = rank[current] + 1
    skip_to(len(names))
    return dict(failed=0, validation=POLICY, tool=tool, cap=cap, reads=len(names),
                intervals=totals['intervals'], zero_count_intervals=totals['zero'],
                reported_locations_checked=totals['locations'], locations_in_bounds=totals['locations'],
                ms_positions_checked=totals['ms'], emitted_counts=totals['intervals'],
                complete_set_intervals=totals['complete'], complete_set_locations=totals['complete_locations'],
                over_cap_intervals=totals['over'], over_cap_locations=totals['over_locations'],
                all_coordinates_checked=False, occurrence_reference_checked=False,
                coordinate_checks='contig bounds; sequence identity by cross-tool agreement at count <= cap, deferred above',
                count_checks='cross-tool agreement')


def read_fasta(path):
    """Read names and sequences in file order; every sequence must be ACGT only."""
    records, name, chunks = {}, None, []
    with Path(path).open() as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    records[name] = ''.join(chunks)
                name, chunks = line[1:].split()[0], []
                if name in records:
                    raise ValueError('Duplicate read identifier')
            else:
                if name is None or line.strip('ACGT'):
                    raise ValueError('Workload reads must be ACGT only')
                chunks.append(line)
    if name is not None:
        records[name] = ''.join(chunks)
    if not records or not all(records.values()):
        raise ValueError('Empty workload or read')
    if list(records) != [f'q{i}' for i in range(len(records))]:
        raise ValueError('Workload read identifiers must be q0, q1, ...')
    return records


def load_truth(truth, offsets, queries):
    bounds = array.array('Q')
    bounds.frombytes(Path(offsets).read_bytes())
    if sys.byteorder != 'little':
        bounds.byteswap()
    names = list(queries)
    if bounds[0] != 0 or len(bounds) < len(names) + 1:
        raise ValueError('Workload offsets differ from reads')
    out = {}
    with Path(truth).open('rb') as stream:
        for i, name in enumerate(names):
            length = bounds[i + 1] - bounds[i]
            if length != len(queries[name]):
                raise ValueError('Read length differs from workload offsets')
            values = array.array('I')
            values.frombytes(stream.read(length * 4))
            if sys.byteorder != 'little':
                values.byteswap()
            if len(values) != length:
                raise ValueError('Truncated matching-statistics truth')
            out[name] = values
        if stream.read(1):
            raise ValueError('Truth longer than the workload')
    return out
