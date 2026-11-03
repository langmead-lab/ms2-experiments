#!/usr/bin/env python3
"""ropebwt3, MONI, PHONI and AugPHONI indexes of one chr22 reference (data/prep/chr22_panel.py output).

These are the build commands of ms-bakery bench/scaling/run.py, which built the frozen chr22-h10 and
chr22-h100 indexes on Sept 23 2026 (records in results/ladder/CELL/TOOL/build-*.command.json on Rockfish):
  ropebwt3   ropebwt3 build -t2 -d -o index.fmd reference.fa   (forward segments; -d adds the reverse strand)
  moni       moni build -r reference.fa -o index -w 10 -p 100 -k -f   (reference.fa is reference-both.fa)
  phoni      phoni build -r reference.txt -w 10 -p 100 -k   (reference.txt is reference-both.txt)
  aug_phoni  aug build -r reference.txt -w 10 -p 100 -k, then pfp_lcp -f reference.txt -w 10, then
             build_aug reference.txt -t compressed
MONI, PHONI and AugPHONI's build drivers are Python scripts in their CMake build directories, run here with
--python. Each tool gets OUT/TOOL/ with its index files, one command.json, stdout and stderr per step, and
build-complete.json. The tool builds themselves (commits, patches and flags) are in TOOLS.md.

Usage: python3 data/prep/legacy_indexes.py REFERENCE_DIR OUT [--tools ropebwt3 moni phoni aug_phoni]
           --ropebwt3 BIN --moni-build DIR --phoni-build DIR --aug-build DIR [--python PYTHON]
"""
import argparse
from pathlib import Path
import shutil
import sys

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from ms2exp.run import common  # noqa: E402

TOOLS = ('ropebwt3', 'moni', 'phoni', 'aug_phoni')
TIMEOUT = 14400


def commands(tool, src, folder, args):
    """The build commands of one tool; copies the reference into folder first, as the frozen builds did."""
    ref, fa = folder/'reference.txt', folder/'reference.fa'
    if tool == 'ropebwt3':
        return [[args.ropebwt3, 'build', '-t2', '-d', '-o', folder/'index.fmd', src/'reference.fa']]
    if tool == 'moni':
        shutil.copy2(src/'reference-both.fa', fa)
        return [[args.python, Path(args.moni_build)/'moni', 'build', '-r', fa, '-o', folder/'index', '-w', '10', '-p', '100',
                 '-k', '-f']]
    shutil.copyfile(src/'reference-both.txt', ref)
    if tool == 'phoni':
        return [[args.python, Path(args.phoni_build)/'phoni', 'build', '-r', ref, '-w', '10', '-p', '100', '-k']]
    base = Path(args.aug_build)
    return [[args.python, base/'aug', 'build', '-r', ref, '-w', '10', '-p', '100', '-k'],
            [base/'_deps/pfp_thresholds-build/test/src/pfp_lcp', '-f', ref, '-w', '10'],
            [base/'test/src/build_aug', ref, '-t', 'compressed']]


def build(tool, src, out, args):
    folder = Path(out)/tool
    folder.mkdir(parents=True, exist_ok=True)
    done = folder/'build-complete.json'
    if done.exists():
        return
    results = [common.execute(c, folder, f'build-{i}', timeout=TIMEOUT)
               for i, c in enumerate(commands(tool, Path(src), folder, args))]
    common.save(done, dict(tool=tool, reference=str(src), commands=results))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('reference', help='a chr22-hN directory of chr22_panel.py')
    p.add_argument('out')
    p.add_argument('--tools', nargs='+', default=list(TOOLS), choices=TOOLS)
    p.add_argument('--ropebwt3')
    p.add_argument('--moni-build')
    p.add_argument('--phoni-build')
    p.add_argument('--aug-build')
    p.add_argument('--python', default=sys.executable)
    args = p.parse_args(argv)
    needed = {'ropebwt3': 'ropebwt3', 'moni': 'moni_build', 'phoni': 'phoni_build', 'aug_phoni': 'aug_build'}
    for t in args.tools:
        if not getattr(args, needed[t]):
            p.error(f'{t} needs --{needed[t].replace("_", "-")}')
    for t in args.tools:
        build(t, args.reference, args.out, args)
        print(t, 'built', flush=True)


if __name__ == '__main__':
    main()
