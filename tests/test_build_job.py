"""The build job's upload of the build records its steps write."""
import json
from pathlib import Path

from ms2exp.build import job
from ms2exp.lib import Store


def test_build_job_uploads_build_records(tmp_path):
    out_root = tmp_path/'outputs'
    (out_root/'tools'/'ms2-abc1234').mkdir(parents=True)
    (out_root/'tools'/'ms2-abc1234'/'ms2-bench').write_text('binary')
    (out_root/'tools'/'build-abc1234.json').write_text('{"compiler": "g++"}')
    builder = job.Builder.__new__(job.Builder)
    builder.results = Store(str(tmp_path/'results'))
    builder.prefix = 'm/build-tools'
    names = builder.publish_provenance(out_root)
    assert names == ['m/build-tools/tools/build-abc1234.json']
    assert json.loads(Path(builder.results.local_path(names[0])).read_text()) == {'compiler': 'g++'}
