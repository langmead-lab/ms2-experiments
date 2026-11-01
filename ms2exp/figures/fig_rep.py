#!/usr/bin/env python3
"""Fig Rep: the cost of each matching statistics reposition in ms2-psi, ms2-phi and ms2-dual.

PLAN.md, October 27: ms2-dual runs the psi and phi LCE walks of a reposition side by side
and stops when either decides it, so its cost per reposition should stay near twice the
smaller of the psi and phi costs. This module measures that per reposition with ms-replog,
a counting build of Orbit (branch reposition-hist, compiled with TMS_REP_LOG), which logs one
64-byte record per reposition (TmsRepRecord in tms_query.hpp). Every mode chooses the same
candidates, so the three modes' logs line up record for record.

Steps:
- simulate (dev node): a batch of PacBio-like reads from ACGT-only source records, with
  template lengths drawn from the real PacBio workload as chr22_sim.py does. The error
  method is one of
    exact: error-free substrings, so repositions come only from differences between the
      query genome and the indexed genomes;
    inject: independent per-base errors at a given rate, split into substitutions,
      insertions and deletions in the proportions pbsim3 3.0.5 gives single-pass
      QSHMM-RSII reads (INJECT_MIX); used below pbsim3's accuracy ceiling;
    pbsim: pbsim3 single pass, QSHMM-RSII, --accuracy-mean 1 - rate, which pbsim3 reads
      to two decimals. pbsim3 ignores accuracies above about 0.966 (a dev probe on October 27 gave the same reads for 0.99,
      0.999 and 0.9999), so this method serves error rates of about 3.4 percent and up.
- run (dev node or Slurm): ms-replog query in modes psi, phi and dual on one batch, one
  thread, interleave 32, as the Fig Ops counting runs; index load time and query time are
  recorded separately.
- summarize: checks that the logs line up and that their sums equal the stats line's
  totals, then writes summary.json with the box statistics of steps per reposition.

Steps per reposition are psi steps plus phi steps in the mode's own run. Repositions with a
match length of 0 take no walk in any mode; they are counted, then left out of the boxes.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import random
import socket
import subprocess
import time

from ms2exp.sim import chr22_sim as sim

MODES = ('psi', 'phi', 'dual')
# Substitution, insertion and deletion shares of pbsim3 single-pass QSHMM-RSII errors at
# --accuracy-mean 0.95 (rates 0.002866, 0.026612, 0.018942 on the dev probe).
INJECT_MIX = (0.0596, 0.5534, 0.3870)
PERCENTILES = (5, 25, 50, 75, 95)


def numpy_dtype():
    import numpy as np
    return np.dtype([('pattern', '<u4'), ('at', '<u4'), ('cap', '<u8'), ('lce', '<u8'), ('dist', '<u8'),
                     ('interval', '<u8'), ('scan', '<u8'), ('psi', '<u8'), ('phi', '<u8')])


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


# Simulation.

def inject(seq, rate, rng, counts):
    """seq with independent errors at rate per base, in INJECT_MIX proportions."""
    if rate <= 0:
        return seq
    out, i, n = [], 0, len(seq)
    while True:
        # Bases until the next error are geometric with parameter rate.
        gap = int(rng.expovariate(-math.log1p(-rate))) if rate < 1 else 0
        if i + gap >= n:
            out.append(seq[i:])
            break
        out.append(seq[i:i + gap])
        i += gap
        u = rng.random()
        if u < INJECT_MIX[0]:
            out.append(rng.choice([b for b in 'ACGT' if b != seq[i]]))
            i += 1
            counts['sub'] += 1
        elif u < INJECT_MIX[0] + INJECT_MIX[1]:
            out.append(rng.choice('ACGT'))
            counts['ins'] += 1
        else:
            i += 1
            counts['del'] += 1
    return ''.join(out)


def simulate(args):
    dest = Path(args.destination).resolve()
    dest.mkdir(parents=True, exist_ok=False)
    source, segments = sim.load_source(args.source, args.pair_reverse_complements, args.source_note)
    real, lengths = sim.real_lengths('pacbio')
    rng = random.Random(args.seed)
    made = sim.templates(segments, lengths, args.bases, rng)
    details = dict(method=args.method, rate=args.rate)
    if args.method == 'exact':
        reads = [(s, o) for _, s, o in made]
    elif args.method == 'inject':
        counts = {'sub': 0, 'ins': 0, 'del': 0}
        reads = [(inject(s, args.rate, rng, counts), o) for _, s, o in made]
        tb = sum(len(s) for _, s, _ in made)
        details.update(mix=dict(zip(('sub', 'ins', 'del'), INJECT_MIX)), counts=counts, template_bases=tb,
                       errors_per_template_base=sum(counts.values()) / tb)
    else:
        work = dest/'work'
        work.mkdir()
        with (work/'templates.fa').open('w') as f:
            for name, seq, _ in made:
                f.write(f'>{name}\n{seq}\n')
        pbsim = sim.binary('pbsim')
        model = sim.MODEL['qshmm-rsii']
        cmd = sim.run_tool([pbsim['path'], '--strategy', 'templ', '--template', work/'templates.fa', '--method', 'qshmm',
                            '--qshmm', model, '--accuracy-mean', f'{1 - args.rate:.2f}', '--pass-num', 1,
                            '--seed', args.seed % (1 << 31), '--prefix', work/'pbsim'], work, 'pbsim')
        origin = {name: o for name, _, o in made}
        template_of = sim.parse_templ_maf(work/'pbsim.maf.gz')
        reads = [(seq, dict(origin[template_of[r]], template=template_of[r]))
                 for r, seq, _ in sim.fastq(work/'pbsim.fq.gz')]
        details.update(tool=pbsim, command=cmd, model=dict(path=str(model), sha256=sim.sha256(model)),
                       accuracy=sim.pbsim_stats(work/'pbsim.stderr'))
    kept = [(s, o) for s, o in reads if s and not set(s) - set('ACGT')]
    with (dest/'reads.fa').open('w') as f:
        for i, (s, _) in enumerate(kept):
            f.write(f'>r{i}\n{s}\n')
    with (dest/'origins.tsv').open('w') as f:
        f.write('read\tsegment\tbegin\tend\tstrand\n')
        for i, (_, o) in enumerate(kept):
            f.write(f"r{i}\t{o['segment']}\t{o['begin']}\t{o['end']}\t{o['strand']}\n")
    manifest = dict(format='fig-rep-reads-v1', created=time.time(), host=socket.gethostname(), seed=args.seed,
                    source=source, lengths_from=real, target_bases=args.bases, details=details,
                    reads=len(kept), dropped=len(reads) - len(kept), bases=sum(len(s) for s, _ in kept),
                    files={'reads.fa': dict(sha256=sha256(dest/'reads.fa'))})
    sim.save_new(dest/'manifest.json', manifest)
    print(json.dumps(dict(reads=manifest['reads'], bases=manifest['bases'], dropped=manifest['dropped'],
                          details={k: v for k, v in details.items() if k not in ('tool', 'command', 'model')})))


# Runs.

def stats_line(stderr):
    rows = []
    for line in stderr.splitlines():
        # Lines named for a mode (stats: dual ...) repeat psi_steps/base; the tms line has no name.
        if line.startswith('stats: ') and 'psi_steps/base=' in line and '=' in line.split()[1]:
            rows.append({k: float(v) for k, _, v in (x.partition('=') for x in line.split()[1:])})
    if len(rows) != 1:
        raise ValueError(f'expected one tms stats line, found {len(rows)}')
    return rows[0]


def fields(stderr, prefix):
    lines = [l for l in stderr.splitlines() if l.startswith(prefix)]
    out = {}
    for item in lines[-1].split()[1:] if lines else []:
        k, eq, v = item.partition('=')
        if eq:
            try:
                out[k] = float(v)
            except ValueError:
                out[k] = v
    return out


def run(args):
    dest = Path(args.destination)
    dest.mkdir(parents=True, exist_ok=False)
    reads = Path(args.reads)
    record = dict(format='fig-rep-run-v1', host=socket.gethostname(), started=time.time(), reads=str(reads),
                  reads_sha256=sha256(reads), binary=args.binary, binary_sha256=sha256(args.binary),
                  interleave=args.interleave, modes={})
    for mode in MODES:
        index = getattr(args, 'index_' + mode)
        log = dest/f'{mode}.replog'
        out = dest/f'{mode}.ms'
        argv = [args.binary, 'query', index, str(reads), '-o', str(out), '--threads', '1',
                '--interleave', str(args.interleave), '--mode', mode]
        started = time.time()
        res = subprocess.run(argv, capture_output=True, text=True, env=dict(os.environ, TMS_REP_LOG=str(log)))
        (dest/f'{mode}.stderr').write_text(res.stderr)
        if res.returncode:
            raise RuntimeError(f'{mode} failed: ' + res.stderr[-2000:])
        q = fields(res.stderr, 'query: ')
        record['modes'][mode] = dict(argv=argv, index=index, wall_seconds=time.time() - started,
                                     index_load_s=q.get('index_load_s'), query_s=q.get('query_s'),
                                     query_ns_per_base=q.get('query_ns_per_base'), total_s=q.get('total_s'),
                                     stats=stats_line(res.stderr), output_sha256=sha256(out))
        out.unlink()
    shas = {m: r['output_sha256'] for m, r in record['modes'].items()}
    if len(set(shas.values())) != 1:
        raise ValueError(f'MS outputs differ across modes: {shas}')
    record['outputs_identical'] = True
    record['finished'] = time.time()
    sim.save_new(dest/'run.json', record)
    summarize_dir(dest)


# Summaries.

def box(values):
    import numpy as np
    q = np.percentile(values, PERCENTILES)
    return dict(zip((f'p{p}' for p in PERCENTILES), map(float, q)), mean=float(values.mean()), n=int(len(values)))


def summarize_dir(dest):
    import numpy as np
    dest = Path(dest)
    dt = numpy_dtype()
    record = json.loads((dest/'run.json').read_text())
    logs = {m: np.fromfile(dest/f'{m}.replog', dtype=dt) for m in MODES}
    # The decision fields must agree in every mode.  interval and scan count LF intervals,
    # which differ between index layouts (the phi layout splits LF differently), so they are
    # only reported.
    key = ('pattern', 'at', 'cap', 'lce', 'dist')
    ref = logs['psi']
    for m in MODES:
        for k in key:
            if len(logs[m]) != len(ref) or not np.array_equal(logs[m][k], ref[k]):
                raise ValueError(f'{m} log differs from psi in {k}')
    same_lf = {m: bool(np.array_equal(logs[m]['interval'], ref['interval']) and
                       np.array_equal(logs[m]['scan'], ref['scan'])) for m in MODES}
    checks = {}
    for m in MODES:
        st = record['modes'][m]['stats']
        bases, reps = st['bases'], len(logs[m])
        got = dict(repositions=reps, psi=int(logs[m]['psi'].sum()), phi=int(logs[m]['phi'].sum()))
        want = dict(repositions=st['repositions/base'] * bases, psi=st['psi_steps/base'] * bases,
                    phi=st['phi_steps/base'] * bases)
        # The stats line prints six significant digits.
        for k in got:
            if abs(got[k] - want[k]) > 1e-5 * max(want[k], 1) + 1:
                raise ValueError(f'{m} {k}: log sum {got[k]} against stats line {want[k]:.0f}')
        checks[m] = dict(got, stats_line=want)
    bases = int(record['modes']['psi']['stats']['bases'])
    walked = ref['cap'] > 0
    cost = {m: (logs[m]['psi'] + logs[m]['phi'])[walked].astype(np.float64) for m in MODES}
    lo = np.minimum(cost['psi'], cost['phi'])
    summary = dict(format='fig-rep-summary-v1', run=str(dest), bases=bases, repositions=int(len(ref)),
                   repositions_per_base=len(ref) / bases, cap0=int((~walked).sum()),
                   aligned='every record matches psi in pattern, at, cap, lce and dist',
                   interval_and_scan_match_psi=same_lf,
                   sums=checks, outputs_identical=record['outputs_identical'],
                   times={m: {k: record['modes'][m][k] for k in ('index_load_s', 'query_s', 'query_ns_per_base')}
                          for m in MODES},
                   steps_per_reposition={m: box(cost[m]) for m in MODES},
                   min_psi_phi=box(lo), ratio_dual_to_min=box(cost['dual'] / lo),
                   dual_le_2min_frac=float((cost['dual'] <= 2 * lo).mean()),
                   dual_le_min_plus1_frac=float((cost['dual'] <= lo + 1).mean()),
                   psi_lt_phi_frac=float((cost['psi'] < cost['phi']).mean()),
                   cap=box(ref['cap'][walked].astype(np.float64)), dist=box(ref['dist'][walked].astype(np.float64)))
    (dest/'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({k: summary[k] for k in ('bases', 'repositions', 'repositions_per_base', 'cap0')}))
    for m in MODES:
        print(m, json.dumps(summary['steps_per_reposition'][m]))
    print('ratio', json.dumps(summary['ratio_dual_to_min']), 'le2min', summary['dual_le_2min_frac'])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('simulate')
    s.add_argument('--source', required=True)
    s.add_argument('--pair-reverse-complements', action='store_true')
    s.add_argument('--source-note')
    s.add_argument('--destination', required=True)
    s.add_argument('--method', required=True, choices=('exact', 'inject', 'pbsim'))
    s.add_argument('--rate', type=float, default=0.0, help='errors per template base')
    s.add_argument('--bases', type=int, default=5_000_000)
    s.add_argument('--seed', type=int, required=True)
    r = sub.add_parser('run')
    r.add_argument('--binary', required=True)
    r.add_argument('--reads', required=True)
    r.add_argument('--destination', required=True)
    for m in MODES:
        r.add_argument('--index-' + m, required=True)
    r.add_argument('--interleave', type=int, default=32)
    u = sub.add_parser('summarize')
    u.add_argument('run_dirs', nargs='+')
    args = ap.parse_args()
    if args.cmd == 'simulate':
        simulate(args)
    elif args.cmd == 'run':
        run(args)
    else:
        for d in args.run_dirs:
            summarize_dir(d)


if __name__ == '__main__':
    main()
