"""Validate default-length SMEM intervals against MS truth and legacy search counts."""
import array
from itertools import groupby
from pathlib import Path
import sys

def groups(path):
    def rows():
        previous = -1
        with Path(path).open() as f:
            for line in f:
                fields = line.split()
                if not fields:
                    continue
                if len(fields) != 4:
                    raise ValueError('Unexpected locating fields in MS comparator')
                name, start, end, count = fields
                i = int(name[1:])
                if name != f'q{i}' or i < previous:
                    raise ValueError('Unexpected or out-of-order read identifier')
                previous = i
                yield i, int(start), int(end), int(count)
    for i, entries in groupby(rows(), key=lambda row: row[0]):
        matches = {}
        for _, start, end, count in entries:
            if start < 0 or end <= start or count <= 0 or (start, end) in matches:
                raise ValueError('Invalid or duplicate SMEM')
            matches[start, end] = count
        yield i, matches


def check_default(output, truth, offsets, baseline=False, count_reference=None, min_length=19):
    n = 1 if baseline else len(offsets) - 1
    observed = iter(groups(output)); current = next(observed, None)
    reference = iter(groups(count_reference)) if count_reference else None
    expected_counts = next(reference, None) if reference is not None else None
    intervals = 0
    with Path(truth).open('rb') as f:
        for i in range(n):
            length = offsets[i + 1] - offsets[i]
            ms = array.array('I'); ms.frombytes(f.read(length * 4))
            if sys.byteorder != 'little': ms.byteswap()
            if len(ms) != length: raise ValueError('Truncated MS truth')
            expected = {(j, j + x) for j, x in enumerate(ms)
                        if x >= min_length and (j == 0 or x >= ms[j - 1])}
            matches = current[1] if current and current[0] == i else {}
            if set(matches) != expected:
                raise ValueError(f'Missing/extra SMEM intervals for q{i}')
            if reference is not None:
                counts = expected_counts[1] if expected_counts and expected_counts[0] == i else {}
                if matches != counts: raise ValueError(f'SMEM occurrence counts differ for q{i}')
                if expected_counts and expected_counts[0] == i: expected_counts = next(reference, None)
            intervals += len(matches)
            if current and current[0] == i: current = next(observed, None)
    if current is not None or expected_counts is not None:
        raise ValueError('Unexpected extra read records')
    return dict(failed=0, reads=n, intervals=intervals, min_length=min_length,
                validation='all long-SMEM intervals from exact MS truth; counts from original MEM search' if reference is not None else 'all long-SMEM intervals from exact MS truth',
                occurrence_counts_checked=reference is not None, exact_ms_equivalence=False)
