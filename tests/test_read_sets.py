"""data/read-sets.json and the code that writes it."""
import io
import json
from pathlib import Path

from ms2exp import read_sets

REPO = Path(__file__).resolve().parents[1]


def test_fasta_stats_and_length_summary():
    sha, size, lengths = read_sets.fasta_stats(io.BytesIO(b'>a\nACGT\nAC\n>b\nA\n>c\nACGTACGT\n'))
    assert size == 28 and lengths == [6, 1, 8]
    assert read_sets.length_summary(lengths) == dict(min=1, median=6, n50=8, max=8)
    assert read_sets.length_summary([1, 2, 3, 4])['median'] == 2.5


def test_every_freeze_read_set_is_recorded_with_its_catalog_checksum():
    sets = json.loads((REPO/'data/read-sets.json').read_text())['sets']
    inputs = json.loads((REPO/'results/freeze-2885ee8-00e665de/manifest.json').read_text())['inputs']
    names = {n[len('reads/'):-len('/reads.fa')] for n in inputs if n.startswith('reads/') and n.endswith('/reads.fa')}
    assert set(sets) == names
    for name, record in sets.items():
        assert record['reads_fa'] == dict(sha256=inputs[f'reads/{name}/reads.fa']['sha256'],
                                          bytes=inputs[f'reads/{name}/reads.fa']['bytes'])
        assert record == dict(read_sets.describe(name.split('/')[1]), **{k: record[k] for k in
                                                                        ('reads_fa', 'reads', 'bases', 'length')})


def test_read_set_counts_agree_with_the_group_summaries():
    sets = json.loads((REPO/'data/read-sets.json').read_text())['sets']
    by_name = {name.split('/')[1]: record for name, record in sets.items()}
    for path in (REPO/'figures/summaries').glob('*/summary.json'):
        w = json.loads(path.read_text())['config']['workload']
        record = by_name.get(w['read_set'])
        if record is None or w.get('reads') is None:
            continue
        assert (record['reads'], record['bases']) == (w['reads'], w['bases']), path.parent.name
