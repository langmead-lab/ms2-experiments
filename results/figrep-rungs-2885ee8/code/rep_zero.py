import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str("/scratch16/blangme2/ms-bakery/bench/campaign-20261012"))
from fig_rep import MODES, numpy_dtype
dt = numpy_dtype()
for run in sys.argv[1:]:
    run = Path(run)
    walked = np.fromfile(run/'psi.replog', dtype=dt)['cap'] > 0
    row = []
    for m in MODES:
        log = np.fromfile(run/f'{m}.replog', dtype=dt)
        c = (log['psi'].astype(np.int64) + log['phi'])[walked]
        row.append(f'{m}: zero={np.mean(c == 0):.3f} mean={c.mean():.2f} p50={np.median(c):.0f}')
        del log, c
    print(f'{run.name:22s}', ' | '.join(row), flush=True)
