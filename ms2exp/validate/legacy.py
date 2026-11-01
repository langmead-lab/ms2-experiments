"""Check every MONI, PHONI or AugPHONI matching-statistics position."""
import array
from pathlib import Path
import sys


def check_lengths(output, truth, offsets, baseline=False):
    """Accept ordered numeric or q-prefixed labels and exact per-read lengths."""
    boundaries = list(offsets)
    if (len(boundaries) < 2 or boundaries[0] != 0 or
            any(a >= b for a, b in zip(boundaries, boundaries[1:]))):
        raise ValueError('Invalid truth boundaries')
    if Path(truth).stat().st_size != boundaries[-1] * 4:
        raise ValueError('Truth size differs from boundaries')
    expected_reads = 1 if baseline else len(boundaries) - 1
    reads = 0
    positions = 0
    current = None
    expected = None
    within = 0
    with Path(truth).open('rb') as reference, Path(output).open() as source:
        for line in source:
            text = line.strip()
            if not text:
                continue
            if text.startswith('>'):
                if current is not None:
                    if within != len(expected):
                        raise ValueError('Incomplete matching-statistics record')
                    reads += 1
                label = text[1:].split()
                if reads >= expected_reads or len(label) != 1 or label[0] not in (str(reads), 'q'+str(reads)):
                    raise ValueError('Unexpected read label or extra record')
                current = reads
                length = boundaries[reads+1] - boundaries[reads]
                expected = array.array('I', reference.read(length*4))
                if sys.byteorder != 'little':
                    expected.byteswap()
                if len(expected) != length:
                    raise ValueError('Truncated truth')
                within = 0
                continue
            if current is None:
                raise ValueError('Values before a read label')
            for token in text.split():
                value = int(token)
                if within >= len(expected):
                    raise ValueError('Extra matching-statistics position')
                if value < 0 or value > len(expected)-within or value != expected[within]:
                    raise ValueError('Matching-statistics mismatch at q%d position %d' % (current, within))
                within += 1
                positions += 1
    if current is not None:
        if within != len(expected):
            raise ValueError('Incomplete matching-statistics record')
        reads += 1
    if reads != expected_reads or positions != boundaries[expected_reads]:
        raise ValueError('Incomplete matching-statistics output')
    return dict(reads=reads, positions=positions, failed=0,
                every_position_checked=True, validation='exact-ms-lengths')
