"""Every number the paper quotes agrees with the committed results, or is a recorded paper edit."""
import copy
import json
from pathlib import Path

from ms2exp.figures import paper_numbers as pn

REPO = Path(__file__).resolve().parents[1]
SPEC = json.loads((REPO/'figures/paper-values.json').read_text())
VALUES, TABLE, NOTES = pn.compute(REPO/'figures', REPO/'results', REPO/'data/read-sets.json')


def statuses(spec):
    return {e['id']: status for e, _, status in pn.check(VALUES, spec)}


def test_strict_check_passes_on_the_committed_results():
    s = statuses(SPEC)
    failed = sorted(k for k, v in s.items() if v == 'FAIL')
    assert not failed, f'the paper and the data disagree on {failed}'
    edits = sorted(e['id'] for e in SPEC['numbers'] if e.get('expected_change'))
    assert all(s[k] in ('edit', 'unchecked') for k in edits)


def test_a_paper_value_that_drifts_from_the_data_fails():
    spec = copy.deepcopy(SPEC)
    entry = next(e for e in spec['numbers'] if e['id'] == 'ops.rec_hprc2_genomes')
    entry['paper'] = '0.00056'
    assert statuses(spec)['ops.rec_hprc2_genomes'] == 'FAIL'


def test_a_recorded_edit_fails_once_the_paper_agrees_or_the_data_move():
    spec = copy.deepcopy(SPEC)
    entry = next(e for e in spec['numbers'] if e['id'] == 'load.short_load_over_query')
    entry['expected_change'] = {'to': entry['paper'], 'edit': 'a recorded edit the paper already carries'}
    assert statuses(spec)['load.short_load_over_query'] == 'FAIL'
    entry['paper'], entry['expected_change']['to'] = '9--23', '9--21'
    assert statuses(spec)['load.short_load_over_query'] == 'FAIL'


def test_rounding_follows_the_paper_text():
    assert pn.agrees(pn.Value('range', [7.74, 13.6]), '7.7--14')
    assert not pn.agrees(pn.Value('range', [7.74, 13.6]), '7.7--13')
    assert pn.agrees(pn.Value('value', 0.95002), '1.0') and pn.agrees(pn.Value('value', 22.498), '22')
    assert pn.agrees(pn.Value('within', [6.07]), '6') and not pn.agrees(pn.Value('within', [7.78]), '6')
    assert pn.agrees(pn.Value('atleast', [22.93]), '22') and not pn.agrees(pn.Value('atleast', [22.93]), '23')


def test_thp_comparison_uses_the_cells_timed_in_both_runs():
    D = pn.Data(REPO/'figures', REPO/'results', REPO/'data/read-sets.json')
    _, _, common, timed, copied = pn.thp_comparison(D)
    assert (len(common), len(timed), len(copied)) == (279, 230, 49)


def test_committed_outputs_are_current():
    rows = pn.check(VALUES, SPEC)
    assert (REPO/'figures/draft/paper-numbers.md').read_text() == pn.report(rows, SPEC, NOTES)
    assert (REPO/'figures/tex/tab-load.tex').read_text() == pn.load_table_tex(TABLE)
    assert r'ropebwt3 & matching statistics & 15 & 4.2 (4.2--4.3) & 3.4 \\' in pn.load_table_tex(TABLE)
