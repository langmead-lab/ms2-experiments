"""Reference and index preparation scripts in data/prep, on tiny fixtures with a stand-in agc."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT/'data/prep'/(name + '.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# Prints the requested interval of a fixed contig in two FASTA lines, as agc getctg does.
FAKE_AGC = '''#!{python}
import sys
contig = 'ACGTNNacgtAAAACCCCGGGGTTTT'
region = sys.argv[-1].split(':')[1]
begin, end = map(int, region.split('-'))
seq = contig[begin:end + 1]
print('>x'); print(seq[:5]); print(seq[5:])
'''

HEAD = 'rank\thaplotype\tcontig\tbegin\tend\torientation\tcontig_md5\tagp_sha256\n'


class Panel(unittest.TestCase):
    def setUp(self):
        self.mod = load('chr22_panel')
        self.dir = tempfile.TemporaryDirectory()
        self.d = Path(self.dir.name)
        (self.d/'order.txt').write_text('100000_CHM13.pri\n200001_HG002.pat\n200002_HG002.mat\n')
        (self.d/'sets.txt').write_text('agc_set\tpansn_prefix\n100000_CHM13.pri\tCHM13#0\n200001_HG002.pat\tHG002#1\n'
                                      '200002_HG002.mat\tHG002#2\n')

    def tearDown(self):
        self.dir.cleanup()

    def test_extract(self):
        panel = self.d/'panel.tsv'
        panel.write_text(HEAD + '1\t200001_HG002.pat\tHG002#1#chr22\t0\t10\t+\tm\ta\n'
                         '2\t200002_HG002.mat\tHG002#2#chr22\t10\t18\t+\tm\tb\n')
        agc = self.d/'agc'
        agc.write_text(FAKE_AGC.format(python=sys.executable))
        agc.chmod(0o755)
        from ms2exp.prep import lists
        check = lists.check_panel
        try:
            lists.check_panel = lambda p=None: check(panel, self.d/'order.txt', self.d/'sets.txt', count=2)
            m = self.mod.extract(2, 'archive.agc', str(agc), self.d/'out')
        finally:
            lists.check_panel = check
        out = self.d/'out/chr22-h2'
        self.assertEqual((out/'reference.fa').read_text(),
                         '>h200001_c0_s0\nACGT\n>h200001_c0_s1\nACGT\n>h200002_c1_s0\nAAAACCCC\n')
        self.assertEqual((out/'reference-both.txt').read_text(), 'ACGTNACGTNAAAACCCCNACGTNACGTNGGGGTTTT')
        self.assertEqual(m['forward_acgt_bases'], 16)
        self.assertEqual(m['sources'][0]['segments'][1], dict(record='h200001_c0_s1', begin=6, end=10))
        self.assertEqual(m['reference_sha256'], hashlib.sha256((out/'reference.fa').read_bytes()).hexdigest())

    def test_agp(self):
        agp = self.d/'HG002.1/ragtag_output/ragtag.scaffold.agp'
        agp.parent.mkdir(parents=True)
        agp.write_text('## agp-version 2.1\n'
                       'CHM13#0#chr21_RagTag\t1\t5\t1\tW\tHG002#1#c9\t1\t5\t+\n'
                       'CHM13#0#chr22_RagTag\t1\t10\t1\tW\tHG002#1#chr22\t1\t10\t+\n'
                       'CHM13#0#chr22_RagTag\t11\t110\t2\tU\t100\tscaffold\tyes\talign_genus\n')
        self.assertEqual(self.mod.agp_rows(agp), [('HG002#1#chr22', 0, 10, '+')])
        self.assertEqual(self.mod.agp_directory('200001_HG002.pat'), 'HG002.1')


class Legacy(unittest.TestCase):
    def test_commands(self):
        mod = load('legacy_indexes')
        with tempfile.TemporaryDirectory() as d:
            src, out = Path(d)/'ref', Path(d)/'out'
            src.mkdir()
            for name in ('reference.fa', 'reference-both.fa', 'reference-both.txt'):
                (src/name).write_text(name)
            args = argparse.Namespace(ropebwt3='rb3', moni_build='M', phoni_build='P', aug_build='A', python='py')
            folder = out/'aug_phoni'
            folder.mkdir(parents=True)
            cmds = mod.commands('aug_phoni', src, folder, args)
            self.assertEqual([c[1] if c[0] == 'py' else c[0] for c in cmds],
                             [Path('A/aug'), Path('A/_deps/pfp_thresholds-build/test/src/pfp_lcp'),
                              Path('A/test/src/build_aug')])
            self.assertEqual((folder/'reference.txt').read_text(), 'reference-both.txt')
            folder = out/'moni'
            folder.mkdir()
            self.assertEqual(mod.commands('moni', src, folder, args)[0][-8:],
                             ['-o', folder/'index', '-w', '10', '-p', '100', '-k', '-f'])
            self.assertEqual((folder/'reference.fa').read_text(), 'reference-both.fa')
            self.assertEqual(mod.commands('ropebwt3', src, out, args)[0][:5],
                             ['rb3', 'build', '-t2', '-d', '-o'])


class Shell(unittest.TestCase):
    def test_syntax(self):
        for script in sorted((ROOT/'data/prep').glob('*.sh')):
            subprocess.run(['bash', '-n', script], check=True)


class Lock(unittest.TestCase):
    def test_prep_tools(self):
        lock = json.loads((ROOT/'tools.lock.json').read_text())
        prep = lock['prep_tools']
        for name, entry in prep.items():
            if isinstance(entry, dict) and entry.get('binary_sha256'):
                self.assertRegex(entry['binary_sha256'], '^[0-9a-f]{64}$', name)
        self.assertEqual(prep['ropebwt3_truth']['binary_sha256'], lock['ropebwt3']['binary_sha256'])


if __name__ == '__main__':
    unittest.main()
