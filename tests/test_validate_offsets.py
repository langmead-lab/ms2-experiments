"""Cross-tool offset policy: parsing, canonical files, comparison, and the workload readers."""
import array
import gzip
import hashlib
from pathlib import Path
import tempfile
import unittest

from ms2exp.validate.offsets import (CanonicalSink, check_smems, compare_files, load_truth, ms_from_smems,
                                     read_canonical, read_contigs, read_fasta)

CONTIGS = {'chr1': 20, 'chr2': 10}
QUERIES = {'q0': 'ACGTAC', 'q1': 'GGA'}
TRUTH = {'q0': [4, 3, 4, 3, 2, 1], 'q1': [3, 2, 1]}
# ms2 query --report smem-all output: one line per SMEM.
ORBIT_CAP1 = 'q0\t0\t4\t3\t1\tchr1:+:5\nq0\t2\t6\t1\t1\tchr2:-:2\nq1\t0\t3\t1\t1\tchr1:+:0\n'
ORBIT_CAP100 = 'q0\t0\t4\t3\t3\tchr1:-:10\tchr2:+:1\tchr1:+:5\nq0\t2\t6\t1\t1\tchr2:-:2\nq1\t0\t3\t1\t1\tchr1:+:0\n'
ROPE_CAP100 = 'q0\t0\t4\t3\t3\tchr1:+:5\tchr1:-:10\tchr2:+:1\nq0\t2\t6\t1\t1\tchr2:-:2\nq1\t0\t3\t1\t1\tchr1:+:0\n'
ROPE_CAP1 = 'q0\t0\t4\t3\t1\tchr2:+:1\nq0\t2\t6\t1\t1\tchr2:-:2\nq1\t0\t3\t1\t1\tchr1:+:0\n'


def check(text, cap=100, tool='orbit', path=None, queries=QUERIES, truth=TRUTH):
    with tempfile.TemporaryDirectory() as tmp:
        output = Path(tmp)/'out';output.write_text(text)
        sink = CanonicalSink(cap, path)
        evidence = check_smems(output, queries, truth, CONTIGS, cap, sink, tool)
        return evidence, sink.close()


def canonical_file(root, text, cap, name, tool='orbit'):
    evidence, summary = check(text, cap, tool, root/name)
    return root/name, summary


class ParserTests(unittest.TestCase):
    def test_orbit_lines_checked_and_counted(self):
        evidence, summary = check(ORBIT_CAP100)
        self.assertEqual((evidence['intervals'], evidence['reported_locations_checked']), (3, 5))
        self.assertEqual((evidence['complete_set_intervals'], evidence['over_cap_intervals']), (3, 0))
        self.assertEqual(evidence['ms_positions_checked'], 9)
        self.assertFalse(evidence['occurrence_reference_checked'])
        evidence, summary = check(ORBIT_CAP1, cap=1)
        self.assertEqual((evidence['over_cap_intervals'], evidence['over_cap_locations']), (1, 1))
        self.assertEqual(summary['over_cap_lines'], 1)

    def test_ms_reconstruction(self):
        self.assertEqual(ms_from_smems([(0, 4), (2, 6)], 6), TRUTH['q0'])
        self.assertEqual(ms_from_smems([], 3), [0, 0, 0])

    def test_structural_failures(self):
        first, second, third = ORBIT_CAP1.splitlines(keepends=True)
        bad = {'text positions without SEQ_IDS': (ORBIT_CAP1.replace('chr1:+:5', '1234'), 1),
               'outside contig': (ORBIT_CAP1.replace('chr1:+:5', 'chr1:+:17'), 1),
               'unknown contig': (ORBIT_CAP1.replace('chr1:+:5', 'chrX:+:5'), 1),
               'missing SMEM': (first + third, 1),
               'extra interval': (ORBIT_CAP1.replace('q0\t2\t6', 'q0\t1\t6'), 1),
               'n_pos': (ORBIT_CAP100.replace('3\t3\tchr1:-:10\tchr2:+:1', '3\t2\tchr2:+:1'), 100),
               'duplicate position': (ORBIT_CAP100.replace('chr2:+:1', 'chr1:+:5'), 100),
               'duplicate interval': (ORBIT_CAP1 + third, 1),
               'read order': (third + first + second, 1)}
        for kind, (text, cap) in bad.items():
            with self.subTest(kind=kind), self.assertRaises(ValueError):
                check(text, cap=cap)

    def test_zero_count_lines_and_silent_reads(self):
        queries, truth = dict(QUERIES, q2='AC'), dict(TRUTH, q2=[0, 1])
        text = ORBIT_CAP1 + 'q2\t0\t1\t0\nq2\t1\t2\t1\t1\tchr1:+:3\n'
        evidence, summary = check(text, 1, queries=queries, truth=truth)
        self.assertEqual(evidence['zero_count_intervals'], 1)
        with self.assertRaises(ValueError):  # A read with SMEMs must have lines.
            check(ORBIT_CAP1, 1, queries=queries, truth=truth)
        queries, truth = dict(queries, q2='A'), dict(truth, q2=[0])
        with self.assertRaises(ValueError):
            check(ORBIT_CAP1 + 'q2\t0\t1\t0\t1\tchr1:+:3\n', 1, queries=queries, truth=truth)

    def test_contig_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'x.len.gz'
            with gzip.open(path, 'wt') as f:f.write('chr1\t20\nchr2\t10\n')
            self.assertEqual(read_contigs(path), CONTIGS)
            with gzip.open(path, 'wt') as f:f.write('chr1\t20\nchr1\t10\n')
            with self.assertRaises(ValueError):read_contigs(path)


class CapTests(unittest.TestCase):
    def test_cap10_is_checked(self):
        evidence, _ = check(ORBIT_CAP100, 10)
        self.assertEqual((evidence['cap'], evidence['complete_set_intervals']), (10, 3))

    def test_nonpositive_cap_refused(self):
        for cap in (0, True, 1.5):
            with self.subTest(cap=cap), self.assertRaises(ValueError):
                check(ORBIT_CAP100, cap)


class LongReadTests(unittest.TestCase):
    def test_long_read_with_many_smems_and_positions(self):
        length, cap = 20000, 100
        smems = [(i, i + 120) for i in range(0, length - 120, 50)] + [(length - 120, length)]
        ms = ms_from_smems(smems, length)
        contigs = {f'c{k}': 10**6 for k in range(cap)}
        lines = ''.join('\t'.join(['q0', str(s), str(e), '472', str(cap)] + [f'c{k}:+:{s}' for k in range(cap)]) + '\n'
                        for s, e in smems)
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)/'out'
            output.write_text(lines)
            sink = CanonicalSink(cap, Path(tmp)/'c.gz')
            evidence = check_smems(output, {'q0': 'A' * length}, {'q0': array.array('I', ms)}, contigs, cap, sink, 'rope')
            summary = sink.close()
        self.assertEqual((evidence['intervals'], evidence['ms_positions_checked']), (len(smems), length))
        self.assertEqual(evidence['over_cap_locations'], cap * len(smems))
        self.assertEqual(summary['lines'], len(smems))


class CanonicalTests(unittest.TestCase):
    def test_sorted_deterministic_and_tool_independent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a, sa = canonical_file(root, ORBIT_CAP100, 100, 'a.gz')
            b, sb = canonical_file(root, ROPE_CAP100, 100, 'b.gz', 'rope')
            self.assertEqual(sa['sha256'], sb['sha256'])
            self.assertEqual(sa['compressed_sha256'], sb['compressed_sha256'])
            lines = gzip.open(a, 'rt').read().splitlines()
            self.assertEqual(lines[0], 'q0\t0\t4\t3\t3\tchr1:+:5\tchr1:-:10\tchr2:+:1')
            records = list(read_canonical(a, sa['sha256']))
            self.assertEqual(records[2], ('q1', 0, 3, 1, (('chr1', '+', 0),)))
            with self.assertRaises(ValueError):list(read_canonical(a, '0'*64))

    def test_over_cap_subset_is_hashed(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, summary = canonical_file(Path(tmp), ORBIT_CAP1, 1, 'c.gz')
            over = ''.join(l for l in gzip.open(path, 'rt') if int(l.split('\t')[3]) > 1)
            self.assertEqual(summary['over_cap_sha256'], hashlib.sha256(over.encode()).hexdigest())
            self.assertEqual((summary['over_cap_lines'], summary['over_cap_locations']), (1, 1))


class ComparisonTests(unittest.TestCase):
    def pair(self, orbit, rope, cap):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a, sa = canonical_file(root, orbit, cap, 'o.gz')
            b, sb = canonical_file(root, rope, cap, 'r.gz', 'rope')
            return compare_files((a, sa['sha256']), (b, sb['sha256']), cap)

    def test_agreement_with_different_order_and_over_cap_choice(self):
        result = self.pair(ORBIT_CAP100, ROPE_CAP100, 100)
        self.assertEqual((result['status'], result['counts']['complete_sets']), ('agree', 3))
        result = self.pair(ORBIT_CAP1, ROPE_CAP1, 1)
        self.assertEqual((result['status'], result['counts']['over_cap_intervals']), ('agree', 1))

    def test_disagreement_on_set_count_and_interval(self):
        result = self.pair(ORBIT_CAP100, ROPE_CAP100.replace('chr2:-:2', 'chr2:-:3'), 100)
        self.assertEqual((result['status'], result['counts']['position_set_differs']), ('disagree', 1))
        self.assertEqual(result['examples'][0]['start'], 2)
        result = self.pair(ORBIT_CAP1, ROPE_CAP1.replace('q0\t0\t4\t3', 'q0\t0\t4\t4'), 1)
        self.assertEqual((result['status'], result['counts']['count_differs']), ('disagree', 1))
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a, sa = canonical_file(root, ORBIT_CAP1, 1, 'o.gz')
            short = root/'s.gz'
            with gzip.GzipFile(short, 'wb', mtime=0) as f:f.write(b'q1\t0\t3\t1\t1\tchr1:+:0\n')
            text = 'q1\t0\t3\t1\t1\tchr1:+:0\n'
            result = compare_files((a, sa['sha256']), (short, hashlib.sha256(text.encode()).hexdigest()), 1)
            self.assertEqual(result['counts']['orbit_only'], 2)


class WorkloadTests(unittest.TestCase):
    """read_fasta and load_truth, which feed check_smems its queries and MS truth."""

    def test_reads_and_truth(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reads = root/'reads.fa';reads.write_text('>q0 first\nACG\nTAC\n\n>q1\nGGA\n')
            queries = read_fasta(reads)
            self.assertEqual(queries, QUERIES)
            truth, offsets = root/'truth.u32', root/'offsets.u64'
            truth.write_bytes(array.array('I', TRUTH['q0'] + TRUTH['q1']).tobytes())
            offsets.write_bytes(array.array('Q', [0, 6, 9]).tobytes())
            self.assertEqual({k: list(v) for k, v in load_truth(truth, offsets, queries).items()}, TRUTH)
            for text in ('>q0\nACGN\n', '>q1\nACG\n', '>q0\nACG\n>q0\nACG\n', '>q0\n', ''):
                with self.subTest(text=text), self.assertRaises(ValueError):
                    reads.write_text(text);read_fasta(reads)
            truth.write_bytes(array.array('I', TRUTH['q0'] + TRUTH['q1'] + [1]).tobytes())
            with self.assertRaises(ValueError):load_truth(truth, offsets, queries)
            truth.write_bytes(array.array('I', TRUTH['q0']).tobytes())
            with self.assertRaises(ValueError):load_truth(truth, offsets, queries)
            truth.write_bytes(array.array('I', TRUTH['q0'] + TRUTH['q1']).tobytes())
            offsets.write_bytes(array.array('Q', [0, 5, 9]).tobytes())
            with self.assertRaises(ValueError):load_truth(truth, offsets, queries)


if __name__ == '__main__':unittest.main()
