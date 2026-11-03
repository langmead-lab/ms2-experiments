"""Simulated chr22 read sets: templates, simulators (faked), read sets and truth."""
import array
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import sys
import tempfile
import unittest
from unittest.mock import patch

from ms2exp.sim import chr22_sim as sim

# Prints one SMEM covering each whole read (count 1), so the MS of position j is length - j.
FAKE_ROPE = '''#!{python}
import sys
name = None
for line in open(sys.argv[-1]):
    line = line.strip()
    if line.startswith('>'):
        name = line[1:]
    elif line:
        print(name, 0, len(line), 1, sep='\\t')
'''
# Writes art.fq and art.aln for two reads; the second holds an N.
FAKE_ART = '''#!{python}
import sys
out = sys.argv[sys.argv.index('-o') + 1]
reads = [('c0_s0-1', 'ACGTA' * 25, 10, '+'), ('c0_s0-2', 'ACGTN' * 25, 20, '-'), ('c0_s0-3', 'TTGCA' * 25, 30, '-')]
with open(out + '.fq', 'w') as fq, open(out + '.aln', 'w') as aln:
    aln.write('##ART_Illumina\\tread_length\\t125\\n##Header End\\n')
    for name, seq, pos, strand in reads:
        fq.write('@' + name + '\\n' + seq + '\\n+\\n' + 'G' * len(seq) + '\\n')
        aln.write('>c0_s0\\t' + name + '\\t' + str(pos) + '\\t' + strand + '\\n' + seq + '\\n' + seq + '\\n')
'''
# pbsim3 multi-pass stand-in: one BAM-named file holding the templates.
FAKE_PBSIM = '''#!{python}
import sys, shutil
a = sys.argv
shutil.copyfile(a[a.index('--template') + 1], a[a.index('--prefix') + 1] + '.bam')
'''
# ccs stand-in: every template except the second becomes a consensus read named S/<zmw>/ccs.
# It refuses to run unless TMPDIR is its working directory, where its temporary files belong.
FAKE_CCS = '''#!{python}
import sys, gzip, os
assert os.path.realpath(os.environ['TMPDIR']) == os.path.realpath(os.getcwd())
names, seqs = [], []
for line in open(sys.argv[-2]):
    line = line.strip()
    if line.startswith('>'):
        names.append(line[1:])
    elif line:
        seqs.append(line)
with gzip.open(sys.argv[-1], 'wt') as out:
    for k, seq in enumerate(seqs, 1):
        if k != 2:
            out.write('@S/%d/ccs np=10 rq=0.999\\n%s\\n+\\n%s\\n' % (k, seq, 'I' * len(seq)))
'''

# pbsim3 template-strategy stand-in: each template t<k> becomes read S_<k+1>, error-free, with its MAF and stats.
FAKE_PBSIM_ONT = '''#!{python}
import sys, gzip
a = sys.argv
prefix = a[a.index('--prefix') + 1]
names, seqs = [], []
for line in open(a[a.index('--template') + 1]):
    line = line.strip()
    if line.startswith('>'):
        names.append(line[1:])
    elif line:
        seqs.append(line)
with gzip.open(prefix + '.fq.gz', 'wt') as fq, gzip.open(prefix + '.maf.gz', 'wt') as maf:
    for k, (name, seq) in enumerate(zip(names, seqs), 1):
        fq.write('@S_%d\\n%s\\n+\\n%s\\n' % (k, seq, '!' * len(seq)))
        maf.write('a\\ns %s 0 %d + %d %s\\ns S_%d 0 %d + %d %s\\n\\n' % (name, len(seq), len(seq), seq, k, len(seq), len(seq), seq))
sys.stderr.write('read num. : %d\\nread accuracy mean (SD) : 0.955301 (0.042975)\\nsubstitution rate. : 0.013539\\n' % len(seqs))
'''


# The held-out haplotype a chr22 manifest names.
SAMPLE = '200245_NA19240.pat'


def tool(folder, name, text):
    path = Path(folder)/name
    path.write_text(text.format(python=sys.executable))
    path.chmod(0o755)
    return path


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class TemplateTests(unittest.TestCase):
    def test_templates_fit_both_strands_and_lengths(self):
        rng = random.Random(4)
        source = ''.join(random.Random(1).choice('ACGT') for _ in range(5000))
        segments = [('a', source[:3000]), ('b', source[3000:3100]), ('c', source[3100:])]
        made = sim.templates(segments, [50, 200, 1500], 40000, rng)
        self.assertGreaterEqual(sum(len(s) for _, s, _ in made), 40000)
        seqs = dict(segments)
        strands = set()
        for name, seq, origin in made:
            self.assertIn(len(seq), (50, 200, 1500))
            part = seqs[origin['segment']][origin['begin']:origin['end']]
            self.assertEqual(seq, part if origin['strand'] == '+' else sim.reverse_complement(part))
            strands.add(origin['strand'])
            self.assertFalse(origin['segment'] == 'b' and len(seq) > 100)
        self.assertEqual(strands, {'+', '-'})
        with self.assertRaises(ValueError):
            sim.templates(segments, [4000], 10, random.Random(1))

    def test_quantiles(self):
        q = sim.quantiles([1, 2, 3, 4, 10])
        self.assertEqual((q['min'], q['median'], q['max'], q['n50'], q['bases']), (1, 3, 10, 10, 20))


class ReadSetTests(unittest.TestCase):
    def finish(self, root, reads, target=10):
        hap = dict(sample=SAMPLE, acgt_bases=100, manifest_sha256='h' * 64)
        with patch.object(sim, 'TARGET', target):
            return sim.finish('ont', reads, root/'set', 5, dict(tools=[]), hap)

    def test_finish_drops_non_acgt_and_stops_at_target(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reads = [('acgt', dict(strand='+')), ('ACNT', dict(strand='-')), ('GGGG', dict(strand='-')),
                     ('TTTTT', dict(strand='+')), ('CCC', dict(strand='+')), ('AAAAAA', dict(strand='-'))]
            m = self.finish(root, reads)
            self.assertEqual(m['simulated']['dropped_non_acgt'], dict(reads=1, bases=4))
            self.assertGreaterEqual(m['bases'], 10)
            offsets = array.array('Q', (root/'set'/'offsets.u64').read_bytes())
            records = list(sim.fasta(root/'set'/'reads.fa'))
            self.assertEqual([n for n, _ in records], [f'q{i}' for i in range(len(records))])
            self.assertEqual(offsets[-1], m['bases'])
            self.assertTrue(all(set(s) <= set('ACGT') for _, s in records))
            # Whole reads only: removing the last would fall below the target.
            self.assertLess(m['bases'] - len(records[-1][1]), 10)
            self.assertEqual((root/'set'/'one.fa').read_text(), f'>q0\n{records[0][1]}\n')
            self.assertEqual(m['files']['reads.fa']['sha256'], digest(root/'set'/'reads.fa'))
            with self.assertRaises(ValueError):
                self.finish(root/'again', reads, target=1000)

    def test_short_simulation_with_fake_art(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            art = tool(root, 'art', FAKE_ART)
            work = root/'work'
            work.mkdir()
            reads, details = sim.simulate_short(root/'acgt.fa', [], work, 9, 1.1, 1000,
                                                art=dict(path=str(art), sha256='x', version='fake'))
            self.assertEqual(len(reads), 3)
            self.assertEqual(reads[1][1], dict(segment='c0_s0', begin=20, end=145, strand='-', source='c0_s0-2'))
            self.assertIn('-ss', details['commands'][0]['argv'])
            self.assertIn('HS25', details['commands'][0]['argv'])

    def test_pacbio_simulation_with_fake_tools(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fake = lambda name, text: dict(path=str(tool(root, name, text)), sha256='x', version='fake')
            work = root/'work'
            work.mkdir()
            source = ''.join(random.Random(2).choice('ACGT') for _ in range(4000))
            with patch.object(sim, 'TARGET', 3000), patch.dict(sim.MODEL, {'qshmm-rsii': root/'art'}):
                (root/'art').write_text('model')
                reads, details = sim.simulate_pacbio([('c0_s0', source)], work, 3, 1.0, [500, 700], 2,
                                                     pbsim=fake('pbsim', FAKE_PBSIM), ccs=fake('ccs', FAKE_CCS),
                                                     samtools=fake('samtools', '#!/bin/sh\n'))
            self.assertEqual(len(reads), details['templates']['count'] - 1)
            self.assertNotIn('accuracy', details)
            self.assertNotIn('t1', [o['template'] for _, o in reads])
            for seq, origin in reads:
                part = source[origin['begin']:origin['end']]
                self.assertEqual(seq, part if origin['strand'] == '+' else sim.reverse_complement(part))


    def test_ccs_accuracy(self):
        a = sim.ccs_accuracy(['I' * 10, '+' * 5 + 'I' * 5])
        self.assertEqual(a['reads'], 2)
        self.assertAlmostEqual(a['median'], 1 - 1e-4)
        self.assertAlmostEqual(a['min'], 1 - (5 * 1e-1 + 5 * 1e-4) / 10)
        self.assertAlmostEqual(a['base_weighted'], 1 - (15 * 1e-4 + 5 * 1e-1) / 20)
        self.assertEqual(sim.ccs_accuracy([]), {})

    def test_pacbio_accuracy_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fake = lambda name, text: dict(path=str(tool(root, name, text)), sha256='x', version='fake')
            work = root/'work'
            work.mkdir()
            source = ''.join(random.Random(3).choice('ACGT') for _ in range(4000))
            with patch.object(sim, 'TARGET', 3000), patch.dict(sim.MODEL, {'qshmm-rsii': root/'model'}):
                (root/'model').write_text('model')
                _, details = sim.simulate_pacbio([('q0', source)], work, 3, 1.0, [500, 700], 2,
                                                 pbsim=fake('pbsim', FAKE_PBSIM), ccs=fake('ccs', FAKE_CCS),
                                                 samtools=fake('samtools', '#!/bin/sh\n'), accuracy=True)
            self.assertAlmostEqual(details['accuracy']['ccs_predicted']['mean_read'], 1 - 1e-4)
            self.assertEqual(details['accuracy']['ccs_predicted']['reads'], details['templates']['count'] - 1)


class OntTests(unittest.TestCase):
    def test_ont_simulation_with_fake_pbsim(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pbsim = dict(path=str(tool(root, 'pbsim', FAKE_PBSIM_ONT)), sha256='x', version='fake')
            work = root/'work'
            work.mkdir()
            source = ''.join(random.Random(5).choice('ACGT') for _ in range(5000))
            (root/'model').write_text('model')
            with patch.object(sim, 'TARGET', 4000), patch.dict(sim.MODEL, {'errhmm-ont-hq': root/'model'}):
                reads, details = sim.simulate_ont([('c0_s0', source)], work, 3, 1.0, [300, 900, 2000], pbsim=pbsim)
            self.assertEqual(len(reads), details['templates']['count'])
            argv = details['commands'][0]['argv']
            self.assertEqual(argv[argv.index('--method') + 1], 'errhmm')
            self.assertNotIn('--sample', argv)
            self.assertEqual(argv[argv.index('--accuracy-mean') + 1], '0.96')
            self.assertEqual(details['accuracy']['read accuracy mean (SD)'], '0.955301 (0.042975)')
            for seq, origin in reads:
                part = source[origin['begin']:origin['end']]
                self.assertEqual(seq, part if origin['strand'] == '+' else sim.reverse_complement(part))
            self.assertEqual({o['strand'] for _, o in reads}, {'+', '-'})


class SourceTests(unittest.TestCase):
    """The same recipe from another source of ACGT-only records (the Fig Ops2 read sets)."""

    def test_load_source_keeps_forward_records_of_pairs(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'src.fa'
            a, b = 'ACGGT' * 20, 'TTGCA' * 30
            path.write_text(f'>q0\n{a}\n>q1\n{sim.reverse_complement(a)}\n>q2\n{b}\n>q3\n{sim.reverse_complement(b)}\n')
            record, segments = sim.load_source(path, pair_reverse_complements=True, note='test')
            self.assertEqual(segments, [('q0', a), ('q2', b)])
            self.assertEqual((record['records'], record['segments'], record['acgt_bases']), (4, ['q0', 'q2'], 250))
            self.assertEqual(record['sha256'], digest(path))
            record, segments = sim.load_source(path)
            self.assertEqual(len(segments), 4)
            path.write_text(f'>q0\n{a}\n>q1\n{a}\n')
            with self.assertRaises(ValueError):
                sim.load_source(path, pair_reverse_complements=True)
            path.write_text('>q0\nACGNT\n')
            with self.assertRaises(ValueError):
                sim.load_source(path)

    def test_ont_from_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pbsim = tool(root, 'pbsim', FAKE_PBSIM_ONT)
            (root/'model').write_text('model')
            real = root/'real'
            real.mkdir()
            (real/'offsets.u64').write_bytes(array.array('Q', [0, 300, 1200, 3200]).tobytes())
            (real/'manifest.json').write_text(json.dumps(dict(
                cohort='real_ont', reads=3, bases=3200, files={'offsets.u64': dict(sha256=digest(real/'offsets.u64'))})))
            source = ''.join(random.Random(6).choice('ACGT') for _ in range(6000))
            path = root/'src.fa'
            path.write_text(f'>q0\n{source}\n>q1\n{sim.reverse_complement(source)}\n')
            loaded = sim.load_source(path, pair_reverse_complements=True, note='held-out test genome')
            with patch.object(sim, 'TARGET', 4000), patch.object(sim, 'P', root/'campaign'), \
                    patch.dict(sim.MODEL, {'errhmm-ont-hq': root/'model'}), \
                    patch.dict(sim.BINARIES, {'pbsim': (pbsim, digest(pbsim), 'fake')}), \
                    patch.dict(sim.REAL, {'ont': real}):
                m = sim.simulate('ont', destination=root/'set', work=root/'work', source=loaded, seed=11,
                                 cohort='fig_ops2_test')
                with self.assertRaises(ValueError):
                    sim.simulate('ont', source=loaded)
                with self.assertRaises(ValueError):
                    sim.simulate('ont', destination=root/'set2')
            self.assertEqual((m['cohort'], m['seed'], m['platform']), ('fig_ops2_test', 11, sim.SETS['ont']['platform']))
            self.assertNotIn('source_haplotype', m)
            self.assertEqual(m['source_segments']['segments'], ['q0'])
            self.assertEqual(m['source_segments']['note'], 'held-out test genome')
            self.assertEqual(m['simulation']['real_workload']['cohort'], 'real_ont')
            self.assertEqual(m['simulation']['accuracy']['read accuracy mean (SD)'], '0.955301 (0.042975)')
            self.assertGreaterEqual(m['bases'], 4000)
            reads = dict(sim.fasta(root/'set'/'reads.fa'))
            rows = [line.split('\t') for line in (root/'set'/'origins.tsv').read_text().splitlines()[1:]]
            for name, length, segment, begin, end, strand, _ in rows:
                part = source[int(begin):int(end)]
                self.assertEqual(segment, 'q0')
                self.assertEqual(reads[name], part if strand == '+' else sim.reverse_complement(part))
            self.assertTrue(list((root/'campaign'/'chr22-sim-jobs').glob('set-local')))

    def test_chr22_manifest_key_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            hap = dict(sample=SAMPLE, acgt_bases=100, manifest_sha256='h' * 64, path='/data/haplotype')
            with patch.object(sim, 'TARGET', 3):
                m = sim.finish('ont', [('ACGT', dict(strand='+'))], Path(tmp)/'set', 5, dict(tools=[]), hap)
            keys = list(m)
            self.assertEqual(keys[keys.index('seed'):keys.index('seed') + 3], ['seed', 'source_haplotype', 'simulation'])
            self.assertEqual(m['cohort'], 'chr22_sim_ont')
            self.assertEqual(m['source_haplotype'], dict(path='/data/haplotype', sample=SAMPLE,
                                                         manifest_sha256='h' * 64, acgt_bases=100))


class TruthTests(unittest.TestCase):
    def setup_reads(self, root, sequences):
        reads = root/'reads'
        reads.mkdir()
        with (reads/'reads.fa').open('w') as f:
            for i, s in enumerate(sequences):
                f.write(f'>q{i}\n{s}\n')
        (reads/'one.fa').write_text(f'>q0\n{sequences[0]}\n')
        offsets = array.array('Q', [0])
        for s in sequences:
            offsets.append(offsets[-1] + len(s))
        (reads/'offsets.u64').write_bytes(offsets.tobytes())
        files = {n: dict(bytes=(reads/n).stat().st_size, sha256=digest(reads/n)) for n in ('reads.fa', 'one.fa', 'offsets.u64')}
        (reads/'manifest.json').write_text(json.dumps(dict(
            status='reads-prepared', set='short', cohort='chr22_sim_short', platform='fake', read_filter=sim.FILTER,
            reads=len(sequences), bases=offsets[-1], one_bases=len(sequences[0]), files=files)))
        return reads

    def test_truth_workload(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reads = self.setup_reads(root, ['ACGTA', 'GGC', 'TTTTTTT'])
            rope = tool(root, 'rope', FAKE_ROPE)
            (root/'index.fmd').write_text('fmd')
            records = lambda reference: (dict(path=str(rope), sha256=digest(rope)),
                                         dict(path=str(root/'index.fmd'), sha256=digest(root/'index.fmd')))
            m = sim.truth('short', 'chr22-h10', reads, root/'workload', root/'work', records=records)
            values = list(array.array('I', (root/'workload'/'truth.u32').read_bytes()))
            self.assertEqual(values, [5, 4, 3, 2, 1, 3, 2, 1, 7, 6, 5, 4, 3, 2, 1])
            self.assertEqual((m['reference'], m['cohort'], m['read_filter']), ('chr22-h10', 'chr22_sim_short', sim.FILTER))
            for name in ('reads.fa', 'one.fa', 'offsets.u64', 'truth.u32'):
                self.assertEqual(m['files'][name]['sha256'], digest(root/'workload'/name))
            argv = json.loads((root/'work'/'truth.command.json').read_text())['argv']
            self.assertEqual(argv[1:4], ['mem', '-t1', '-l1'])
            with self.assertRaises(FileExistsError):
                sim.truth('short', 'chr22-h10', reads, root/'workload', root/'work2', records=records)

    def test_truth_refuses_zero_and_overlong_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'t.u32'
            path.write_bytes(array.array('I', [2, 1, 0]).tobytes())
            with self.assertRaises(ValueError):
                sim.bounded(path, [0, 3])
            path.write_bytes(array.array('I', [4, 2, 1]).tobytes())
            with self.assertRaises(ValueError):
                sim.bounded(path, [0, 3])
            path.write_bytes(array.array('I', [3, 2, 1]).tobytes())
            sim.bounded(path, [0, 3])


class RopeRecordTests(unittest.TestCase):
    def test_rope_records_come_from_the_data_directory(self):
        with patch.dict(os.environ, {'MS2EXP_DATA': '/data/ms2exp'}):
            rope, fmd = sim.rope_records('chr22-h10')
        self.assertEqual(rope['path'], '/data/ms2exp/tools/ropebwt3-cfad8475/ropebwt3')
        self.assertEqual(rope['sha256'], sim.ROPE_SHA256)
        self.assertEqual((fmd['path'], fmd['bytes']), ('/data/ms2exp/ropebwt3/chr22-h10/index.fmd', 104856784))


if __name__ == '__main__':
    unittest.main()
