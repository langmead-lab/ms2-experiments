"""Thread-scaling support: the threads matrices, CPU binding, the cross-thread check, the inputs
catalog from an earlier run, and the 1 Gbase superset read sets."""
import argparse
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from ms2exp import catalog, collect, manifest
from ms2exp.build import d13_prep
from ms2exp.lib import MANIFEST_FORMAT
from ms2exp.run import runner

COMMIT = 'b' * 40
ROPE_SHA256 = 'cfad8475911f9298bf53177bf4a578a6d133f0464f4b9970bef9afb04697d53d'


def matrix(fn):
    b = manifest.Builder({}, True, COMMIT, ROPE_SHA256[:8])
    fn(b)
    return b


def test_threads_matrix_has_seventy_cells_and_an_smt_check_first():
    b = matrix(manifest.threads_matrix)
    assert b.jobs[0]['id'] == 'threads-smt-check'
    smt = b.jobs[0]['cells']
    assert len(smt) == 1 and smt[0]['threads'] == 64 and smt[0]['binding'] == 'smt'
    sweep = [c for j in b.jobs[1:] for c in j['cells']]
    assert len(b.jobs) == 11 and len(sweep) == 70
    assert all(j['platforms'] == ['aws'] for j in b.jobs)
    for c in sweep:
        assert c['hash_output'] and c['figures'] == ['threads']
        assert c['binding'] == ('smt' if c['threads'] == 64 else 'physical')
        assert c['read_set'].endswith('-ms-1gb-acgt-v1')
        if c['tool'] == 'ms2':
            assert c['block_bytes'] == 2097152 and c['mode'] == 'ordinary' and c['interleave'] == 32
        else:
            assert c['tool'] == 'ropebwt3' and c['min_len'] == 1 and 'batch' not in c
    ids = [c['id'] for c in sweep]
    assert len(set(ids)) == 70
    assert 'threads-human_short/phi-ordinary-t16' in ids and 'threads-human_ont_modern/rope-l1-t64' in ids
    # Each job is one configuration and read set at every thread count.
    for job in b.jobs[1:]:
        assert sorted(c['threads'] for c in job['cells']) == list(manifest.THREADS)


def test_pilot_matrix_runs_every_configuration_at_1_and_32_on_the_60mb_sets():
    b = matrix(manifest.threads_pilot_matrix)
    assert [j['id'] for j in b.jobs] == ['threads-pilot'] and b.jobs[0]['platforms'] == ['rockfish']
    cells = b.jobs[0]['cells']
    assert len(cells) == 20
    assert {c['read_set'] for c in cells} == {'human_short-ms-60mb-acgt-v2', 'human_ont_modern-ms-60mb-acgt-v2'}
    assert sorted({c['threads'] for c in cells}) == [1, 32]
    # The 60 Mb sets keep the Fig 5 truth and inputs.
    assert 'reads/human472/human_short-ms-60mb-acgt-v2/truth.u32' in b.inputs


def test_threads_manifest_records_its_protocol(tmp_path, monkeypatch):
    names = sorted(set(matrix(manifest.threads_matrix).inputs))
    inputs = {n: dict(sha256=ROPE_SHA256 if n.endswith('/ropebwt3') else f'{k:064x}', bytes=k + 1)
              for k, n in enumerate(names)}
    (tmp_path/'inputs.json').write_text(json.dumps(dict(inputs=inputs)))
    monkeypatch.setattr(manifest, 'now', lambda: 1.5e9)
    monkeypatch.setattr(sys, 'argv', ['manifest', '--matrix', 'threads', '--inputs', str(tmp_path/'inputs.json'),
                                      '--ms2-repo', 'https://example.org/ms2.git', '--ms2-commit', COMMIT,
                                      '-o', str(tmp_path/'m.json')])
    manifest.main()
    m = json.loads((tmp_path/'m.json').read_text())
    assert m['id'].startswith('threads-' + COMMIT[:7] + '-')
    assert m['protocol']['threads'] == 'per cell' and m['protocol']['numa']['membind'] == 0
    assert m['protocol']['runs'] == 3 and m['protocol']['transparent_hugepages'] == 'always'


def fake_sysfs(root, nodes, smt):
    """nodes: list of core counts per node; hardware thread j of core i is cpu i + j * total cores."""
    total = sum(nodes)
    first = 0
    for n, count in enumerate(nodes):
        cpus = [c + j * total for c in range(first, first + count) for j in range(smt)]
        node = root/'node'/f'node{n}'
        node.mkdir(parents=True)
        (node/'cpulist').write_text(','.join(map(str, sorted(cpus))) + '\n')
        for c in range(first, first + count):
            siblings = [c + j * total for j in range(smt)]
            for cpu in siblings:
                t = root/'cpu'/f'cpu{cpu}'/'topology'
                t.mkdir(parents=True, exist_ok=True)
                (t/'thread_siblings_list').write_text(','.join(map(str, siblings)) + '\n')
        first += count


def test_physical_cores_lists_node_0_first(tmp_path):
    fake_sysfs(tmp_path, [4, 4], 2)
    cores = runner.physical_cores(tmp_path)
    assert cores[:2] == [(0, [0, 8]), (0, [1, 9])]
    assert [n for n, _ in cores] == [0] * 4 + [1] * 4
    assert runner.cpu_list('0-2,5,7-8\n') == [0, 1, 2, 5, 7, 8]


def bare_runner(tmp_path, platform, nodes, smt):
    fake_sysfs(tmp_path/'sys', nodes, smt)
    r = object.__new__(runner.Runner)
    r.platform = platform
    r.args = SimpleNamespace(no_numa=False)
    r.context = {'cores': runner.physical_cores(tmp_path/'sys')}
    r.staged = {'bin': Path('/x/ms2'), 'idx': Path('/x/i'), 'reads': Path('/x/r.fa'), 'fmd': Path('/x/f.fmd')}
    return r


def cell(**kw):
    return dict(dict(tool='ms2', kind='ms', layout='phi', mode='ordinary', interleave=32, block_bytes=2097152,
                     inputs=dict(binary='bin', index='idx', reads='reads', fmd='fmd')), **kw)


def test_binding_on_aws_uses_socket_0_and_both_hardware_threads_at_64(tmp_path):
    r = bare_runner(tmp_path, 'aws', [32, 32], 2)
    assert r.cpus(cell(threads=4, binding='physical')) == [0, 1, 2, 3]
    assert r.cpus(cell(threads=32, binding='physical')) == list(range(32))
    assert r.cpus(cell(threads=64, binding='smt')) == [c for i in range(32) for c in (i, i + 64)]
    with pytest.raises(runner.CellFailed):
        r.cpus(cell(threads=33, binding='physical'))
    argv = r.ms2_argv(cell(threads=8, binding='physical'), Path('/o'))
    assert argv[:3] == ['numactl', '--physcpubind=0,1,2,3,4,5,6,7', '--membind=0']
    assert argv[argv.index('--threads') + 1] == '8' and argv[argv.index('--block-bytes') + 1] == '2097152'
    assert r.env(cell(threads=8))['OMP_NUM_THREADS'] == '8'


def test_binding_on_a_rockfish_pilot_node_spills_onto_node_1(tmp_path):
    r = bare_runner(tmp_path, 'rockfish', [24, 24], 1)
    assert r.cpus(cell(threads=32, binding='physical')) == list(range(32))
    rope = r.rope_argv(cell(tool='ropebwt3', min_len=1, threads=32, binding='physical'))
    assert rope[1].startswith('--physcpubind=0,1,') and '-t32' in rope and '-l1' in rope


def test_cells_without_threads_keep_the_freeze_binding(tmp_path):
    r = bare_runner(tmp_path, 'aws', [2], 1)
    c = cell()
    argv = r.ms2_argv(c, Path('/o'))
    assert argv[:3] == ['numactl', '--cpunodebind=0', '--membind=0']
    assert argv[argv.index('--threads') + 1] == '1' and r.env(c) is None
    assert '-t1' in r.rope_argv(cell(tool='ropebwt3', min_len=1))


def summary(cell_id, threads, digests, tool='ms2', layout='psi', status='checked'):
    return dict(cell=cell_id, status=status, tool=tool, layout=layout, mode='ordinary', min_len=1,
                read_set='human_short-ms-1gb-acgt-v1', threads=threads, output_sha256=digests)


def test_cross_thread_check():
    same = [summary(f'g/psi-ordinary-t{t}', t, ['a']) for t in (1, 2, 64)]
    rope = [summary(f'g/rope-l1-t{t}', t, ['c'], tool='ropebwt3') for t in (1, 2)]
    out = collect.cross_thread(same + rope + [dict(summary('g/x', 1, ['z']), status='failed')])
    assert [(o['config'], o['threads'], o['passes']) for o in out] == [
        ('ms2-psi-ordinary', [1, 2, 64], True), ('ropebwt3-1', [1, 2], True)]
    differ = collect.cross_thread(same + [summary('g/psi-ordinary-t4', 4, ['a', 'b'])])
    assert not differ[0]['passes'] and differ[0]['sha256'] == ['a', 'b']


def test_catalog_from_results_resolves_built_inputs(tmp_path):
    results = tmp_path/'run'
    (results/'build-x').mkdir(parents=True)
    (results/'manifest.json').write_text(json.dumps(dict(format=MANIFEST_FORMAT, id='freeze-x', inputs={
        'reads/a': dict(sha256='1' * 64, bytes=5), 'ms2/i': dict(sha256='built', bytes=None, producer='build-x'),
        'traces/t': dict(sha256='pending', bytes=None)})))
    (results/'build-x'/'outputs.json').write_text(json.dumps(dict(outputs={'ms2/i': dict(sha256='2' * 64, bytes=9)})))
    extra = tmp_path/'extra.json'
    extra.write_text(json.dumps(dict(inputs={'reads/b': dict(sha256='3' * 64, bytes=7)})))
    out = tmp_path/'inputs.json'
    catalog.cmd_from_results(argparse.Namespace(results=str(results), add=[str(extra)], out=str(out)))
    inputs = json.loads(out.read_text())['inputs']
    assert sorted(inputs) == ['ms2/i', 'reads/a', 'reads/b']
    assert inputs['ms2/i']['sha256'] == '2' * 64 and inputs['ms2/i']['bytes'] == 9
    extra.write_text(json.dumps(dict(inputs={'reads/a': dict(sha256='4' * 64, bytes=5)})))
    with pytest.raises(SystemExit, match='differ'):
        catalog.cmd_from_results(argparse.Namespace(results=str(results), add=[str(extra)], out=str(out)))


def fake_pool(root):
    """Six pool reads; r1 and r4 hold an N. Core is r0, r2, r4 (r4 split at its N), so the filtered
    core order is r0, r2 and the rest r3, r5."""
    reads = ['ACGTACGTAC', 'ACGNACGTAC', 'TTTTGGGGCC', 'acgtacgtac', 'AAAANCCCCC', 'GGGGGGGGGG']
    root.mkdir(parents=True)
    (root/'selected.fa').write_text(''.join(f'>r{i}\n{s}\n' for i, s in enumerate(reads)))
    (root/'original-read-names.tsv').write_text(''.join(f'r{i}\t1\torig{i}\t{len(s)}\n' for i, s in enumerate(reads)))
    (root/'core.map.tsv').write_text('q0\tr0\t0\t10\nq1\tr2\t0\t10\nq2\tr4\t0\t4\nq3\tr4\t5\t10\n')


def prefix_set(root, seqs):
    root.mkdir(parents=True)
    (root/'reads.fa').write_text(''.join(f'>q{i}\n{s}\n' for i, s in enumerate(seqs)))
    d13_prep.write_u64(root/'offsets.u64', [0] + [10 * (i + 1) for i in range(len(seqs))])
    (root/'one.fa').write_text(f'>q0\n{seqs[0]}\n')
    (root/'manifest.json').write_text(json.dumps(dict(reads=len(seqs), bases=10 * len(seqs))))


def superset(tmp_path, monkeypatch, name, prefix, limit):
    monkeypatch.setattr(d13_prep, 'POOLS', {'human_short': tmp_path/'pool'})
    monkeypatch.setattr(d13_prep, 'OUT', tmp_path/'out')
    d13_prep.cmd_superset(argparse.Namespace(name=name, cohort='human_short', prefix=str(prefix), limit=limit))
    return tmp_path/'out'/name


def test_superset_puts_the_core_first_and_keeps_the_prefix(tmp_path, monkeypatch):
    fake_pool(tmp_path/'pool')
    prefix_set(tmp_path/'p60', ['ACGTACGTAC'])
    dest = superset(tmp_path, monkeypatch, 's', tmp_path/'p60', 35)
    seqs = [line for line in (dest/'reads.fa').read_text().splitlines() if not line.startswith('>')]
    # Core first (r0, r2), then the rest in pool order (r3 uppercased); r5 would pass 35 bases.
    assert seqs == ['ACGTACGTAC', 'TTTTGGGGCC', 'ACGTACGTAC']
    assert list(d13_prep.read_u64(dest/'offsets.u64')) == [0, 10, 20, 30]
    m = json.loads((dest/'manifest.json').read_text())
    assert (m['reads'], m['bases'], m['core_reads'], m['rest_reads']) == (3, 30, 2, 1)
    assert m['dropped_non_acgt_reads'] == dict(core=1, rest=1) and m['reached_limit']
    assert (dest/'map.tsv').read_text().splitlines()[3] == 'q2\tr3\torig3\trest'


def test_superset_refuses_a_set_that_is_not_its_prefix(tmp_path, monkeypatch):
    fake_pool(tmp_path/'pool')
    prefix_set(tmp_path/'p60', ['TTTTGGGGCC'])
    with pytest.raises(ValueError, match='not a prefix'):
        superset(tmp_path, monkeypatch, 's', tmp_path/'p60', 100)
