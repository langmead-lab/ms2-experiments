"""figures/draft/CAPTIONS.md states the protocol, read sets and names of the committed records."""
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
CAPTIONS = REPO/'figures/draft/CAPTIONS.md'
# Wording from earlier campaigns that the freeze records contradict.
STALE = ['warmup run.', 'timed runs after one warmup', 'NA19240', 'ms2-lcp', 'layout', 'campaign-20261012',
         'bench/campaign', 'Orbit commit', 'HPRC2 full', 'Rockfish']


def test_captions_carry_no_stale_wording():
    text = CAPTIONS.read_text()
    found = [w for w in STALE if w in text]
    assert not found, f'CAPTIONS.md still says {found}'


def test_captions_state_the_manifest_protocol_and_query_genome():
    text = CAPTIONS.read_text()
    p = json.loads((REPO/'results/freeze-2885ee8-00e665de/manifest.json').read_text())['protocol']
    assert (p['runs'], p['warmup_runs'], p['report']) == (3, 0, 'median')
    assert 'Each value is the median of three runs, with no warmup run' in text
    assert f"transparent huge pages {p['transparent_hugepages']}" in text
    assert 'simulated from APR008 haplotype 1' in text


def test_protocol_sentence_follows_the_records_and_refuses_a_mix():
    pytest.importorskip('matplotlib')
    from ms2exp.figures import paper_figures as pf
    summary = json.loads((REPO/'figures/summaries/human_short-ms-acgt60mb/summary.json').read_text())
    assert pf.protocol_sentence([summary]).startswith('Each value is the median of three runs, with no warmup run')
    other = json.loads(json.dumps(summary))
    other['config']['protocol'].update(runs=2, warmup_runs=1)
    assert pf.protocol_sentence([other]).startswith('Each value is the median of two runs after one warmup run')
    with pytest.raises(ValueError):
        pf.protocol_sentence([summary, other])
