"""The committed input lists and their checks (ms2exp.prep.lists), on the real lists and tiny fixtures."""
import gzip
from pathlib import Path
import tempfile
import unittest

from ms2exp.prep import lists


def write_len(path, names):
    with gzip.open(path, 'wt') as f:
        for name in names:
            f.write(f'{name}\t100\n')


class Haplotypes(unittest.TestCase):
    def test_pansn(self):
        self.assertEqual(lists.pansn('100000_CHM13.pri'), 'CHM13#0')
        self.assertEqual(lists.pansn('200001_HG002.pat'), 'HG002#1')
        self.assertEqual(lists.pansn('200002_HG002.mat'), 'HG002#2')
        self.assertEqual(lists.pansn('210019_APR008.hap1'), 'APR008#1')

    def test_committed_list(self):
        rows = lists.haplotypes()
        self.assertEqual(len(rows), 472)
        self.assertEqual(rows[0]['pansn_prefix'], 'CHM13#0')

    def test_len_file(self):
        with tempfile.TemporaryDirectory() as d:
            sets = Path(d)/'sets.txt'
            sets.write_text('agc_set\tpansn_prefix\n100000_CHM13.pri\tCHM13#0\n200001_HG002.pat\tHG002#1\n'
                            '200002_HG002.mat\tHG002#2\n')
            good = Path(d)/'good.len.gz'
            write_len(good, ['CHM13#0#chr1', 'CHM13#0#chr2', 'HG002#1#c1', 'HG002#2#c1', 'HG002#2#c2'])
            self.assertEqual(lists.check_haplotypes(good, sets, expected=3), 3)
            missing = Path(d)/'missing.len.gz'
            write_len(missing, ['CHM13#0#chr1', 'HG002#1#c1'])
            with self.assertRaises(ValueError):
                lists.check_haplotypes(missing, sets, expected=3)
            reordered = Path(d)/'reordered.len.gz'
            write_len(reordered, ['HG002#1#c1', 'CHM13#0#chr1', 'HG002#2#c1'])
            with self.assertRaises(ValueError):
                lists.check_haplotypes(reordered, sets, expected=3)
            with self.assertRaises(ValueError):
                lists.check_haplotypes(good, sets, expected=472)


class Panel(unittest.TestCase):
    def test_committed_panel(self):
        rows = lists.check_panel()
        self.assertEqual(len(rows), 242)
        self.assertEqual(sum(int(r['end']) - int(r['begin']) for r in rows), 5068343876)
        self.assertEqual(rows[0]['haplotype'], '200001_HG002.pat')
        self.assertEqual(lists.donor_order()[99], '200100_HG01975.mat')

    def test_wrong_rank(self):
        with tempfile.TemporaryDirectory() as d:
            order = Path(d)/'order.txt'
            order.write_text('100000_CHM13.pri\n200001_HG002.pat\n200002_HG002.mat\n')
            sets = Path(d)/'sets.txt'
            sets.write_text('agc_set\tpansn_prefix\n100000_CHM13.pri\tCHM13#0\n200001_HG002.pat\tHG002#1\n'
                            '200002_HG002.mat\tHG002#2\n')
            head = 'rank\thaplotype\tcontig\tbegin\tend\torientation\tcontig_md5\tagp_sha256\n'
            panel = Path(d)/'panel.tsv'
            panel.write_text(head + '1\t200001_HG002.pat\tHG002#1#chr22\t0\t10\t+\tx\ty\n'
                             '2\t200002_HG002.mat\tHG002#2#chr22\t0\t10\t+\tx\ty\n')
            self.assertEqual(len(lists.check_panel(panel, order, sets, count=2)), 2)
            panel.write_text(head + '1\t200002_HG002.mat\tHG002#2#chr22\t0\t10\t+\tx\ty\n'
                             '2\t200001_HG002.pat\tHG002#1#chr22\t0\t10\t+\tx\ty\n')
            with self.assertRaises(ValueError):
                lists.check_panel(panel, order, sets, count=2)


if __name__ == '__main__':
    unittest.main()
