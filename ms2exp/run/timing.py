"""Index load and query phases of one timed run, read from the tool's own phase timer, and the
counting build's operation counts.

Policy load-query-separate-v2: query time runs from the end of index loading through the last
output written, and timing gates apply to it. Load time is reported as its own quantity. Whole-run
wall and CPU time stay recorded as cross-checks.

ms2 query prints one summary line, starting `query:` and carrying total_s=, with index_load_s,
ftab_s, query_s, read_s, parse_s, format_s and total_s in seconds. With an ftab it also prints an
`ftab:` line whose load_s repeats ftab_s. index_load_s covers reading the index and, in oracle
mode, the replay trace (and the counting table smem-all builds); ftab_s covers reading the ftab.
total_s starts after both loads and after the read and output files are opened, and ends after
the last block is read, parsed, queried, formatted and written and the output is flushed. Load is
index_load_s + ftab_s; query is total_s, which like ropebwt3's interval includes reading and
parsing the queries. The parser also accepts the equivalent `batch:` line of older builds.

ropebwt3 stamps each log line with seconds since process start. Load is the last
`rb3_fmi_load_all` "loaded" stamp; query is the last `worker_pipeline` "processed" stamp minus
that. The stamps have millisecond resolution, so a one-read run often shows 0.

MONI, PHONI and AugPHONI print no phase timer; the runner times them by the paired one-read run
(runner.py).
"""
import math
import re

# Modes whose stats lines carry no reposition rate and are merged by mode name.
SUMMARY_LAYOUTS = ('phisum', 'dualsum')
POLICY = 'load-query-separate-v2'
ORBIT_QUERY = 'total_s'
ROPE_QUERY = 'last processed stamp - last loaded stamp'
_ROPE = re.compile(r'^\[M::(\w+)::([0-9.]+)\*[0-9.]+\] (.*)$')


def _seconds(value):
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ValueError('Invalid phase time')
    return value


def orbit(text):
    batch = [line for line in text.splitlines() if line.startswith('batch: ')
             or (line.startswith('query: ') and ' total_s=' in line)]
    if len(batch) != 1:
        raise ValueError('Orbit phase timer line missing or repeated')
    fields = dict(item.split('=', 1) for item in batch[0].split()[1:] if '=' in item)
    parts = {key: _seconds(fields[key]) for key in
             ('index_load_s', 'ftab_s', 'query_s', 'read_s', 'parse_s', 'format_s', 'total_s')}
    ftab = [line for line in text.splitlines() if line.startswith('ftab: ')]
    if ftab:
        loaded = dict(item.split('=', 1) for item in ftab[-1].split()[1:] if '=' in item)
        if not math.isclose(_seconds(loaded['load_s']), parts['ftab_s'], rel_tol=1e-6, abs_tol=1e-9):
            raise ValueError('Orbit ftab load times differ')
    return dict(source='orbit-' + batch[0].split(':', 1)[0] + '-line', query_definition=ORBIT_QUERY,
                load_seconds=parts['index_load_s'] + parts['ftab_s'],
                query_seconds=parts['total_s'], components=parts)


def ropebwt3(text):
    loaded, processed, sequences = [], [], 0
    for line in text.splitlines():
        match = _ROPE.match(line)
        if not match:
            continue
        stage, stamp, message = match.groups()
        if stage == 'rb3_fmi_load_all' and message.startswith('loaded '):
            loaded.append(_seconds(stamp))
        elif stage == 'worker_pipeline' and message.startswith('processed '):
            processed.append(_seconds(stamp))
            sequences += int(message.split()[1])
    if not loaded or not processed or processed[-1] < loaded[-1]:
        raise ValueError('ropebwt3 phase timestamps missing')
    return dict(source='ropebwt3-log-stamps', query_definition=ROPE_QUERY,
                load_seconds=loaded[-1], query_seconds=processed[-1] - loaded[-1],
                components=dict(last_loaded_s=loaded[-1], last_processed_s=processed[-1],
                                processed_sequences=sequences))
PARSERS = {'orbit': orbit, 'ropebwt3': ropebwt3}


def parse(tool, text):
    return PARSERS[tool](text)


def parse_stats(stderr, layout):
    """The counting build's `stats:` line for the layout's index kind, with operations per base.

With both MS_STATS and TMS_STATS, a query prints a line for each kind; the other kind's
line has zero counts (nan or inf ratios). An Orbit-ms line has walk_up/rep; a tms line
has psi_steps/base.
"""
    want = 'walk_up/rep' if layout == 'ms' else 'psi_steps/base'
    rows, summary = [], {}
    for line in stderr.splitlines():
        if not line.startswith('stats: '):
            continue
        fields = {}
        for item in line.split()[1:]:
            key, eq, value = item.partition('=')
            if eq:
                fields[key] = float(value)
        # Lines named for a mode (stats: dual ..., stats: phisum ..., stats: dualsum ...) have no
        # reposition rate, and a dual or dualsum line repeats psi_steps/base. Each mode's lines
        # are merged into one record.
        kind = line.split()[1] if len(line.split()) > 1 else ''
        if kind in SUMMARY_LAYOUTS or (kind and '=' not in kind):
            summary.setdefault(kind, {}).update(fields)
        elif want in fields:
            rows.append(fields)
    if len(rows) != 1:
        raise ValueError(f'Expected one stats line with {want}, found {len(rows)}')
    f = rows[0]
    rep = f['repositions/base']
    if not math.isfinite(rep):
        raise ValueError('Reposition rate is not finite')
    finite = lambda v: v if math.isfinite(v) else 0.0
    if layout == 'ms':
        walk = rep * (finite(f['walk_up/rep']) + finite(f['walk_down/rep'])) if rep else 0.0
        ops = dict(lf=1.0, psi=0.0, phi=0.0, walk=walk, total=1.0 + walk, repositions=rep)
    else:
        psi, phi = finite(f['psi_steps/base']), finite(f['phi_steps/base'])
        ops = dict(lf=1.0, psi=psi, phi=phi, walk=0.0, total=1.0 + psi + phi, repositions=rep,
                   scan_rows_per_base=rep * finite(f.get('scan_rows/rep', 0.0)),
                   walk_visits_per_base=finite(f.get('walk_visits/base', 0.0)))
    ops['lf_ff_per_step'] = finite(f.get('lf_ff/step', 0.0))
    out = dict(line=f, ops_per_base=ops, lf_basis='one LF step per base, as in the TeraMS Figure 5 reproduction')
    if summary:
        out['summary_lines'] = summary
    return out
