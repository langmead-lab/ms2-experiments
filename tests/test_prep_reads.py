"""Read-set preparation (ms2exp.prep and data/download/download_skin.py) on tiny fixtures."""
import array
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from ms2exp.prep import acgt_workload, pieces, pool
from ms2exp.sim import chr22_sim
from ms2exp.sim.chr22_sim import sha256

ROOT = Path(__file__).resolve().parents[1]


def u64(path):
    a = array.array('Q')
    a.frombytes(Path(path).read_bytes())
    return list(a)


def write_fastq(path, reads):
    with gzip.open(path, 'wt') as f:
        for name, seq in reads:
            f.write(f'@{name}\n{seq}\n+\n{"I" * len(seq)}\n')


class Pool(unittest.TestCase):
    def test_selected_and_split(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            write_fastq(d/'s1.fq.gz', [('a/1', 'acgtNNac'), ('b/1', 'GGGG')])
            write_fastq(d/'s2.fq.gz', [('a/2', 'TTTT'), ('b/2', 'CCNN')])
            pool.check_mates(d/'s1.fq.gz', d/'s2.fq.gz')
            self.assertEqual(pool.write_selected([d/'s1.fq.gz', d/'s2.fq.gz'], d), 4)
            self.assertEqual((d/'selected.fa').read_text(), '>r0\nACGTNNAC\n>r1\nGGGG\n>r2\nTTTT\n>r3\nCCNN\n')
            self.assertEqual((d/'original-read-names.tsv').read_text().splitlines()[2], 'r2\t2\ta/2\t4')
            stats = pool.split(d/'selected.fa', d, 'core')
            self.assertEqual((stats['reads'], stats['bases'], stats['ambiguous_bases_removed']), (5, 16, 4))
            self.assertEqual(u64(d/'core.offsets.u64'), [0, 4, 6, 10, 14, 16])
            self.assertEqual((d/'core.map.tsv').read_text().splitlines()[1], 'q1\tr0\t6\t8')

    def test_mates_out_of_order(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            write_fastq(d/'s1.fq.gz', [('a/1', 'A'), ('b/1', 'C')])
            write_fastq(d/'s2.fq.gz', [('b/2', 'A'), ('a/2', 'C')])
            with self.assertRaises(ValueError):
                pool.check_mates(d/'s1.fq.gz', d/'s2.fq.gz')


class AcgtWorkload(unittest.TestCase):
    def test_select(self):
        records = [('q0', 'ACGT'), ('q1', 'ACNT'), ('q2', 'GG'), ('q3', 'TTTT')]
        bounds = [0, 4, 8, 10, 14]
        kept, skipped = acgt_workload.select(records, bounds, lambda i: i != 2, 6)
        self.assertEqual([i for i, _ in kept], [0, 3])
        self.assertEqual((skipped['non_acgt'], skipped['not_whole']), (1, 1))
        with self.assertRaises(ValueError):
            acgt_workload.select(records, bounds, lambda i: True, 100)


class AcgtSplit(unittest.TestCase):
    def test_split(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            with gzip.open(d/'in.fa.gz', 'wt') as f:
                f.write('>s1 x\nacgtNN\nAC\n>s2\nNNNN\n>s3\nGT\n')
            subprocess.run([sys.executable, '-m', 'ms2exp.prep.acgt_split', d/'in.fa.gz', d/'out.fa', d/'off.u64',
                            d/'stats.json'], check=True, cwd=ROOT, capture_output=True)
            self.assertEqual((d/'out.fa').read_text(), '>s1_0_4\nACGT\n>s1_6_8\nAC\n>s3_0_2\nGT\n')
            self.assertEqual(u64(d/'off.u64'), [0, 4, 6, 8])
            stats = json.loads((d/'stats.json').read_text())
            self.assertEqual((stats['input_records_fully_dropped'], stats['non_acgt_bases_dropped']), (1, 6))


@unittest.skipUnless(importlib.util.find_spec('numpy'), 'pieces needs numpy')
class Pieces(unittest.TestCase):
    def test_workload_and_pieces(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d/'src.fa').write_text('>g1\nACGTACGTAC\nNNGG\n>g2\nTTTTT\n')
            w = d/'w'
            m = pieces.make_workload('mtb152', 'tiny', d/'src.fa', w, block=4)
            self.assertEqual((m['reads'], m['bases'], m['dropped_non_acgt']), (3, 17, 2))
            self.assertEqual((w/'reads.fa').read_text(), '>q0\nACGTACGTAC\n>q1\nGG\n>q2\nTTTTT\n')
            # Truth as if every whole record occurs in the index: MS[i] = length - i.
            offsets = u64(w/'offsets.u64')
            values = [n - j for a, b in zip(offsets, offsets[1:]) for n in [b - a] for j in range(n)]
            (w/'truth.u32').write_bytes(array.array('I', values).tobytes())
            files = {n: dict(sha256=sha256(w/n)) for n in pieces.WORKLOAD_FILES}
            (w/'manifest.json').write_text(json.dumps(dict(status='prepared', read_filter=pieces.FILTER,
                                                           reference='mtb152', bases=17, files=files)))
            r = pieces.pieces(w, 4, d/'p')
            self.assertEqual(r['reads'], 6)
            self.assertEqual(u64(d/'p'/'offsets.u64'), [0, 4, 8, 10, 12, 16, 17])
            truth = list(array.array('I', (d/'p'/'truth.u32').read_bytes()))
            self.assertEqual(truth[:10], [4, 3, 2, 1, 4, 3, 2, 1, 2, 1])
            self.assertEqual((d/'p'/'map.tsv').read_text().splitlines()[2], 'q2\tq0\t8\t10')


class Chr22SimCli(unittest.TestCase):
    def test_source_needs_destination(self):
        with self.assertRaises(SystemExit):
            chr22_sim.main(['simulate', '--set', 'pacbio', '--source', 'x.fa'])


class DownloadSkin(unittest.TestCase):
    def test_table(self):
        spec = importlib.util.spec_from_file_location('download_skin', ROOT/'data/download/download_skin.py')
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        files = mod.files(runs=['ERR14066769'])
        self.assertEqual([f['relative_path'] for f in files],
                         ['skin/ERR14066769/ERR14066769_1.fastq.gz', 'skin/ERR14066769/ERR14066769_2.fastq.gz'])
        self.assertEqual(len({f['run'] for f in mod.files()}), 77)
        with self.assertRaises(ValueError):
            mod.files(runs=['ERR0'])


if __name__ == '__main__':
    unittest.main()
