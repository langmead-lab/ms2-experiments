"""The toy manifest builds from a catalog of its own input names, and its id depends only on content."""
import json
import sys

import pytest

from ms2exp import manifest

COMMIT = 'a' * 40
ROPE_SHA256 = 'cfad8475911f9298bf53177bf4a578a6d133f0464f4b9970bef9afb04697d53d'


def toy_names():
    """Every input name the toy matrix asks for, from a Builder with an empty catalog."""
    b = manifest.Builder({}, False, COMMIT, ROPE_SHA256[:8])
    manifest.toy_matrix(b)
    return sorted(set(b.missing))


def catalog(path, names):
    inputs = {}
    for k, name in enumerate(names):
        digest = ROPE_SHA256 if name.endswith('/ropebwt3') else f'{k:064x}'
        inputs[name] = dict(sha256=digest, bytes=k + 1)
    path.write_text(json.dumps(dict(inputs=inputs)))
    return path


def build(tmp_path, monkeypatch, inputs, name, created=1.5e9):
    monkeypatch.setattr(manifest, 'now', lambda: created)
    out = tmp_path/name
    monkeypatch.setattr(sys, 'argv', ['manifest', '--matrix', 'toy', '--inputs', str(inputs),
                                      '--ms2-repo', 'https://example.org/ms2.git', '--ms2-commit', COMMIT,
                                      '-o', str(out)])
    manifest.main()
    return json.loads(out.read_text())


def test_toy_manifest_builds_with_a_stable_id(tmp_path, monkeypatch):
    names = toy_names()
    assert names
    inputs = catalog(tmp_path/'inputs.json', names)
    first = build(tmp_path, monkeypatch, inputs, 'a.json')
    assert first['matrix'] == 'toy' and first['created'] == 1.5e9
    assert first['id'].startswith('toy-' + COMMIT[:7] + '-')
    # Oracle traces may stay pending, since a Rockfish run records them; every other input is cataloged.
    pending = {k for k, v in first['inputs'].items() if v['sha256'] == 'pending'}
    assert pending and all(k.startswith('traces/') and k.endswith('.oracle') for k in pending)
    assert set(first['inputs']) - pending == set(names)
    cells = [c['id'] for j in first['jobs'] for c in j['cells']]
    assert cells and len(cells) == len(set(cells))
    assert build(tmp_path, monkeypatch, inputs, 'b.json')['id'] == first['id']
    assert build(tmp_path, monkeypatch, inputs, 'c.json', created=1.6e9)['id'] != first['id']


def test_missing_input_is_refused(tmp_path, monkeypatch):
    names = toy_names()
    inputs = catalog(tmp_path/'inputs.json', names[1:])
    with pytest.raises(SystemExit, match='missing from the catalog'):
        build(tmp_path, monkeypatch, inputs, 'a.json')
