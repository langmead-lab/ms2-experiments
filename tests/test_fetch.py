"""ms2exp.fetch: sources in order, sha256 checks, duplicates fetched once, and the inputs table."""
import csv
import hashlib
import json

import pytest

from ms2exp import fetch


def write_table(path, files):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fetch.FIELDS, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for name, data in files.items():
            w.writerow(dict(name=name, bytes=len(data), sha256=hashlib.sha256(data).hexdigest(),
                            source='zenodo:' + fetch.DEPOSIT['reads'], file=fetch.deposit_file(name), url='-'))


@pytest.fixture
def setup(tmp_path):
    files = {'reads/a/reads.fa': b'>q0\nACGT\n', 'reads/a/truth.u32': b'\1\0\0\0' * 4,
             'reads/b/reads.fa': b'>q0\nACGT\n'}
    table = tmp_path/'inputs.tsv'
    write_table(table, files)
    mirror = tmp_path/'mirror'
    mirror.mkdir()
    for name, data in files.items():
        (mirror/fetch.deposit_file(name)).write_bytes(data)
    return tmp_path, table, mirror, files


def run(tmp_path, table, *extra):
    fetch.main(['--data', str(tmp_path/'data'), '--table', str(table), *extra])


def test_mirror_fetch_verifies_and_links_duplicates(setup):
    tmp_path, table, mirror, files = setup
    run(tmp_path, table, '--mirror', str(mirror), *files)
    data = tmp_path/'data'
    for name, content in files.items():
        assert (data/name).read_bytes() == content
    assert (data/'reads/a/reads.fa').stat().st_ino == (data/'reads/b/reads.fa').stat().st_ino
    index = json.loads((data/'locations.json').read_text())['locations']
    assert sorted(index) == sorted(files)
    run(tmp_path, table, '--verify-only', *files)


def test_corrupt_source_is_refused(setup):
    tmp_path, table, mirror, files = setup
    (mirror/fetch.deposit_file('reads/a/truth.u32')).write_bytes(b'\2\0\0\0' * 4)
    with pytest.raises(SystemExit, match='sha256'):
        run(tmp_path, table, '--mirror', str(mirror), 'reads/a/truth.u32')
    assert not (tmp_path/'data/reads/a/truth.u32').exists()
    assert not (tmp_path/'data/reads/a/truth.u32.part').exists()


def test_locations_link_and_verify_only(setup):
    tmp_path, table, mirror, files = setup
    lab = tmp_path/'lab'
    lab.mkdir()
    (lab/'r.fa').write_bytes(files['reads/a/reads.fa'])
    locs = tmp_path/'locs.json'
    locs.write_text(json.dumps(dict(locations={'reads/a/reads.fa': str(lab/'r.fa')})))
    with pytest.raises(SystemExit, match='Missing'):
        run(tmp_path, table, '--verify-only', 'reads/a/reads.fa')
    run(tmp_path, table, '--from-locations', str(locs), '--link', 'reads/a/reads.fa')
    dest = tmp_path/'data/reads/a/reads.fa'
    assert dest.is_symlink() and dest.resolve() == (lab/'r.fa').resolve()
    (lab/'r.fa').write_bytes(b'>q0\nACGA\n')
    with pytest.raises(SystemExit, match='differs'):
        run(tmp_path, table, '--verify-only', 'reads/a/reads.fa')


def test_no_public_copy_is_an_error(setup):
    tmp_path, table, mirror, files = setup
    with pytest.raises(SystemExit, match='no public copy'):
        run(tmp_path, table, 'reads/a/reads.fa')


def test_dry_run_reports_inputs_without_a_source(setup, capsys):
    tmp_path, table, mirror, files = setup
    with pytest.raises(SystemExit, match='1 inputs have no public copy'):
        run(tmp_path, table, '--dry-run', 'reads/a/reads.fa', 'reads/b/reads.fa', 'reads/a/truth.u32',
            '--from-locations', str(_locs(tmp_path, files, 'reads/a/truth.u32')))
    out = capsys.readouterr().out
    assert 'NONE    reads/a/reads.fa' in out and 'reads/b/reads.fa' not in out.split('NONE')[-1].split('\n')[0]
    run(tmp_path, table, '--dry-run', '--mirror', str(mirror), *files)


def _locs(tmp_path, files, name):
    lab = tmp_path/'lab-dry'
    lab.mkdir(exist_ok=True)
    (lab/'f').write_bytes(files[name])
    locs = tmp_path/'locs-dry.json'
    locs.write_text(json.dumps(dict(locations={name: str(lab/'f')})))
    return locs


def test_committed_table_covers_the_toy():
    table = fetch.load_table(fetch.TABLE)
    names = [n.strip() for n in fetch.SETS['toy'].read_text().splitlines() if n.strip() and not n.startswith('#')]
    assert len(names) == 34 and all(n in table for n in names)
    for n in names:
        assert table[n]['sha256'] not in ('built', 'pending') and table[n]['file'] != '-'
    fmd = table['ropebwt3/mtb152/index.fmd']
    assert fmd['source'] == 'zenodo:' + fetch.DEPOSIT['precursors']
    assert fmd['url'] == fetch.ZENODO_URL.format(record=fetch.DEPOSIT['precursors'], file=fmd['file'])
    assert all(r['source'] == fetch.UNDECIDED for n, r in table.items() if n.startswith('tools/'))


def test_committed_table_names_the_published_deposits():
    table = fetch.load_table(fetch.TABLE)
    for name, row in table.items():
        assert 'PLACEHOLDER' not in row['source'] + row['url'], name
        if row['source'] in fetch.DEPOSITED:
            record = row['source'].split(':', 1)[1]
            assert row['url'] == fetch.ZENODO_URL.format(record=record, file=row['file']), name
    counts = {record: sum(r['source'] == 'zenodo:' + record for r in table.values())
              for record in fetch.DEPOSIT.values()}
    assert counts == {'23964218': 71, '23964220': 30, '23964222': 18}
