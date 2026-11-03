#!/usr/bin/env python3
"""Matching-statistics truth for the read pools of ms2exp.prep.pool, as made on Rockfish before the freeze.

For each set of a cohort (by default only core, the pool every frozen real-read set is sliced from), this
runs ropebwt3 mem -t THREADS -l1 on the human472 FMD, rebuilds per-position matching statistics from the
SMEMs with ms2exp.validate.ms.check, and writes truth/COHORT/SET.u32 and SET.json (passed, command, and
the set's record from complete.json). The SMEM output is deleted once converted. The frozen runs used
24 threads and the ropebwt3 binary with sha256 cfad8475; their records are
data/evidence/rockfish/results/real-hprc2-1d580b6/truth/COHORT/core.json. Like the pools, this is
provenance and is not rerun by data/prep/run_d13.sh.

Usage: python -m ms2exp.prep.truth_pool COHORT --pools DIR --truth DIR --ropebwt3 PATH --fmd human472.fmd
           [--sets core ...] [--threads 24]
"""
import argparse
import array
import json
from pathlib import Path
import sys

from ms2exp.run import common
from ms2exp.validate.ms import check

ROPE_SHA256 = 'cfad8475911f9298bf53177bf4a578a6d133f0464f4b9970bef9afb04697d53d'
SETS = ('pilot', 'diagnostic', 'offset', 'core', 'scaling')


def read_offsets(path):
    a = array.array('Q')
    a.frombytes(Path(path).read_bytes())
    if sys.byteorder != 'little':
        a.byteswap()
    return a


def truth(cohort, setname, pools, out, rope, fmd, threads=24, timeout=14400):
    inputs = Path(pools)/cohort
    done = Path(out)/cohort/(setname + '.json')
    if done.exists():
        return json.loads(done.read_text())
    offsets = read_offsets(inputs/(setname + '.offsets.u64'))
    folder = Path(out)/cohort/setname
    row = common.execute([rope, 'mem', '-t', threads, '-l1', fmd, inputs/(setname + '.fa')], folder, 'truth', timeout)
    row['check'] = check('ropebwt3', folder, row, Path(out)/cohort/(setname + '.u32'), offsets, write=True)
    Path(row['stdout']).unlink()
    row['raw_output_retained'] = False
    common.save(folder/'truth.command.json', row)
    record = dict(passed=True, command=row,
                  input=json.loads((inputs/'complete.json').read_text())['sets'][setname])
    common.save(done, record)
    return record


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cohort')
    p.add_argument('--pools', required=True, help='the --out directory of ms2exp.prep.pool')
    p.add_argument('--truth', required=True, help='destination; COHORT/SET.u32 and SET.json are written there')
    p.add_argument('--ropebwt3', required=True)
    p.add_argument('--fmd', required=True)
    p.add_argument('--sets', nargs='+', default=['core'], choices=SETS)
    p.add_argument('--threads', type=int, default=24)
    args = p.parse_args(argv)
    if common_sha256(args.ropebwt3) != ROPE_SHA256:
        print('Warning: this ropebwt3 is not the binary that made the frozen truth', file=sys.stderr)
    for s in args.sets:
        r = truth(args.cohort, s, args.pools, args.truth, args.ropebwt3, args.fmd, args.threads)
        print(s, json.dumps(r['command']['check']), flush=True)


def common_sha256(path):
    from ms2exp.sim.chr22_sim import sha256
    return sha256(path)


if __name__ == '__main__':
    main()
