#!/usr/bin/env python3
"""Build or collect the pinned tool binaries of a freeze run into an inputs tree, on Rockfish.

Usage:
  python -m ms2exp.build.tools [--lock tools.lock.json] [--ms2-repo URL --ms2-commit SHA] --root INPUTS_DIR --ropebwt3 PATH
                 [--ropebwt3-sha256 SHA] [--legacy moni=PATH --legacy phoni=PATH --legacy augphoni=PATH]
                 [--subdir src] [--work DIR] [--parts ms2,pinned] [--submodule-repo NAME URL ...]

--parts ms2 builds only what depends on the ms2 commit (the first four files below); the AWS build
job of an all-AWS run uses it. --parts pinned only collects the pinned binaries, which then never
need the commit's source.

Writes, under INPUTS_DIR (logical names as ms2exp.manifest uses them):
  tools/ms2-<c7>/ms2-bench    make ms2-bench at the commit (assertions off)
  tools/ms2-<c7>/ms2-stats    the counting build: ms2-bench plus -DMS_STATS -DTMS_STATS (D7); the
                              Makefile's own ms2-stats target sets only -DMS_STATS, so this
                              build compiles the Makefile's sources itself
  tools/ms2-<c7>/ms2-replog   make ms2-replog: the counting build plus -DTMS_REP_LOG (A1)
  tools/ms_check              campaign-20261012/ms_check.cpp
  tools/ropebwt3-<s8>/ropebwt3  the pinned binary, checked against its sha256
  tools/runtime/libgomp.so.1  the OpenMP runtime the pinned ropebwt3 links, for hosts without it
  tools/<tool>/<binary>       MONI, PHONI and AugPHONI binaries, copied as given
  tools/runtime/libstdc++.so.6, tools/runtime/libgcc_s.so.1
                              the C++ runtime the legacy binaries resolve through ldd (Rockfish's
                              GCC 13), for hosts whose own is older or absent
  tools/build-<c7>.json       provenance: commands, compiler, sha256 of every binary

The ms2 source is fetched by commit from --ms2-repo (a URL, a local repository path, or a git
bundle such as make_source_bundles.sh writes), with its Orbit and TeraTools submodules at the
commits it records. --submodule-repo NAME URL fetches submodule NAME (orbit, teratools) from URL,
which may also be a path or a bundle, in place of the URL in .gitmodules. The C++ binaries link libstdc++ and
libgcc statically so that one set of bytes runs on Rockfish and on the AWS boxes (Amazon Linux
2023 ships an older libstdc++ than Rockfish's GCC 13).
"""
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

from ms2exp.lib import atomic_json, file_sha256

HERE = Path(__file__).resolve().parent

STATIC = '-static-libstdc++ -static-libgcc'
LEGACY_BINARY = {'moni': 'moni-ms', 'phoni': 'phoni', 'augphoni': 'aug_phoni'}
# The C++ runtime the legacy binaries link, which Amazon Linux 2023 may lack; the AWS runs use these copies.
LEGACY_RUNTIME = ('libstdc++.so.6', 'libgcc_s.so.1')


def run(argv, cwd, record, shell=False):
    started = time.time()
    result = subprocess.run(argv, cwd=cwd, capture_output=True, text=True, shell=shell)
    record.append(dict(argv=argv, cwd=str(cwd), returncode=result.returncode, seconds=time.time() - started,
                       log=(result.stdout + result.stderr)[-4000:]))
    if result.returncode:
        raise SystemExit(f'Command failed: {argv}\n{result.stdout}{result.stderr}')
    return result.stdout


def fetch(repo, commit, dest, record, overrides=()):
    """Shallow checkout of one commit and its submodules; GitHub serves any commit reachable from a ref."""
    if (dest/'.git').exists():
        head = run(['git', 'rev-parse', 'HEAD'], dest, record).strip()
        if head == commit:
            submodules(dest, record, overrides)
            return
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    run(['git', 'init', '-q'], dest, record)
    run(['git', 'remote', 'add', 'origin', repo], dest, record)
    run(['git', 'fetch', '-q', '--depth', '1', 'origin', commit], dest, record)
    run(['git', 'checkout', '-q', 'FETCH_HEAD'], dest, record)
    head = run(['git', 'rev-parse', 'HEAD'], dest, record).strip()
    if head != commit:
        raise SystemExit(f'Fetched {head}, not {commit}')
    submodules(dest, record, overrides)


def submodules(dest, record, overrides=()):
    """Orbit (orbit/) and TeraTools (teratools/) at the commits the ms2 commit records."""
    run(['git', 'submodule', 'init'], dest, record)
    for name, url in overrides:
        run(['git', 'config', f'submodule.{name}.url', url], dest, record)
        if (dest/name/'.git').exists():  # a checkout from an earlier attempt keeps its own remote
            run(['git', 'remote', 'set-url', 'origin', url], dest/name, record)
    run(['git', '-c', 'protocol.file.allow=always', 'submodule', 'update', '--depth', '1'], dest, record)


def install(src, dest, provenance, key):
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dest)
    dest.chmod(0o755)
    provenance['binaries'][key] = dict(path=str(dest), sha256=file_sha256(dest))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--lock', help='tools.lock.json: the ms2 repo and commit, the ropebwt3 sha256 and the legacy '
                   'binaries\' sha256, used where the flags below are not given')
    p.add_argument('--ms2-repo')
    p.add_argument('--ms2-commit')
    p.add_argument('--root', required=True, help='inputs tree to write tools/ into')
    p.add_argument('--ropebwt3', help='the pinned ropebwt3 binary (needed with --parts pinned)')
    p.add_argument('--ropebwt3-sha256')
    p.add_argument('--legacy', action='append', default=[], metavar='TOOL=PATH')
    p.add_argument('--subdir', default='src', help='directory of the ms2 Makefile in the repository')
    p.add_argument('--cxx', default='g++')
    p.add_argument('--work', default=None, help='checkout directory (default ROOT/../ms2-src-<c7>)')
    p.add_argument('--parts', default='ms2,pinned', help='comma-separated: ms2, pinned')
    p.add_argument('--submodule-repo', action='append', default=[], nargs=2, metavar=('NAME', 'URL'),
                   help='fetch submodule NAME from URL in place of its .gitmodules URL')
    args = p.parse_args()
    lock = json.loads(Path(args.lock).read_text()) if args.lock else {}
    args.ms2_repo = args.ms2_repo or lock.get('ms2', {}).get('repo')
    args.ms2_commit = args.ms2_commit or lock.get('ms2', {}).get('commit')
    args.ropebwt3_sha256 = (args.ropebwt3_sha256 or lock.get('ropebwt3', {}).get('binary_sha256')
                            or 'cfad8475911f9298bf53177bf4a578a6d133f0464f4b9970bef9afb04697d53d')
    args.legacy_sha256 = {tool: lock[tool]['binary_sha256'] for tool in LEGACY_BINARY if tool in lock}
    if not args.ms2_repo or not args.ms2_commit:
        raise SystemExit('Give --ms2-repo and --ms2-commit, or --lock')
    parts = set(args.parts.split(','))
    if parts - {'ms2', 'pinned'} or not parts:
        raise SystemExit('--parts takes ms2, pinned or both')
    if 'pinned' in parts and not args.ropebwt3:
        raise SystemExit('--parts pinned needs --ropebwt3')
    commit, short = args.ms2_commit, args.ms2_commit[:7]
    if len(commit) != 40:
        raise SystemExit('--ms2-commit must be a full 40-hex commit')
    root = Path(args.root).resolve()
    work = Path(args.work or root.parent/f'ms2-src-{short}').resolve()
    commands = []
    provenance = dict(ms2=dict(repo=args.ms2_repo, commit=commit, subdir=args.subdir), started=time.time(),
                      commands=commands, binaries={}, parts=sorted(parts))
    if 'ms2' in parts:
        build_ms2(args, commit, short, root, work, commands, provenance)
    if 'pinned' in parts:
        collect_pinned(args, root, work, commands, provenance)
    provenance['finished'] = time.time()
    atomic_json(root/f"tools/build-{short}{'' if parts == {'ms2', 'pinned'} else '-' + '-'.join(sorted(parts))}.json",
                provenance)
    print(json.dumps(provenance['binaries'], indent=1))


def build_ms2(args, commit, short, root, work, commands, provenance):
    fetch(args.ms2_repo, commit, work, commands, args.submodule_repo)
    provenance['ms2']['submodules'] = run(['git', 'submodule', 'status'], work, commands).splitlines()
    provenance['ms2']['submodule_repos'] = args.submodule_repo
    source = work/args.subdir
    cxx = f'{args.cxx} {STATIC}'
    provenance['compiler'] = run([args.cxx, '--version'], source, commands).splitlines()[0]
    run(['make', f'CXX={cxx}', f'MS_COMMIT={commit}', 'ms2-bench'], source, commands)
    install(source/'ms2-bench', root/f'tools/ms2-{short}/ms2-bench', provenance, 'ms2-bench')
    # The counting build mirrors the exploratory counting build: the ms2-bench command plus both stats flags.
    makefile = (source/'Makefile').read_text()
    srcs = re.search(r'^MS_SRCS = (.+)$', makefile, re.M).group(1).split()
    stats = source/'ms2-stats-freeze'
    run([args.cxx, *STATIC.split(), '-std=c++17', '-I../orbit/include', '-I.', '-O3', '-pthread',
         f'-DMS_COMMIT="{commit}"', '-DNDEBUG', '-DMS_STATS', '-DTMS_STATS', '-o', str(stats), *srcs, '-lz'],
        source, commands)
    install(stats, root/f'tools/ms2-{short}/ms2-stats', provenance, 'ms2-stats')
    # A1's build: the Makefile target adds -DTMS_REP_LOG to the stats flags (TMS_REP_LOG=path logs repositions).
    run(['make', f'CXX={cxx}', f'MS_COMMIT={commit}', 'ms2-replog'], source, commands)
    install(source/'ms2-replog', root/f'tools/ms2-{short}/ms2-replog', provenance, 'ms2-replog')
    checker = work/'ms_check'
    run([args.cxx, *STATIC.split(), '-O2', '-std=c++17', '-o', str(checker), str(HERE/'ms_check.cpp')],
        work, commands)
    install(checker, root/'tools/ms_check', provenance, 'ms_check')
    provenance['help'] = run([str(root/f'tools/ms2-{short}/ms2-bench'), '--help'], work, commands).splitlines()[0]


def collect_pinned(args, root, work, commands, provenance):
    work.mkdir(parents=True, exist_ok=True)
    if file_sha256(args.ropebwt3) != args.ropebwt3_sha256:
        raise SystemExit('ropebwt3 binary differs from its pinned sha256')
    install(Path(args.ropebwt3), root/f'tools/ropebwt3-{args.ropebwt3_sha256[:8]}/ropebwt3', provenance, 'ropebwt3')
    gomp = [line.split('=>')[1].split()[0] for line in run(['ldd', args.ropebwt3], work, commands).splitlines()
            if 'libgomp' in line and '=>' in line]
    if gomp:
        install(Path(gomp[0]), root/'tools/runtime/libgomp.so.1', provenance, 'libgomp')
    runtime = {}
    for item in args.legacy:
        tool, path = item.split('=', 1)
        if tool in args.legacy_sha256 and file_sha256(path) != args.legacy_sha256[tool]:
            raise SystemExit(f'{tool} binary differs from its pinned sha256 in the lock file')
        install(Path(path), root/f'tools/{tool}/{LEGACY_BINARY[tool]}', provenance, tool)
        for line in run(['ldd', path], work, commands).splitlines():
            soname, _, target = line.partition('=>')
            if soname.strip() in LEGACY_RUNTIME and target.split():
                runtime.setdefault(soname.strip(), set()).add(target.split()[0])
    for soname, paths in sorted(runtime.items()):
        if len(paths) > 1:
            raise SystemExit(f'Legacy binaries resolve {soname} to different files: {sorted(paths)}')
        install(Path(paths.pop()), root/f'tools/runtime/{soname}', provenance, soname)
    if args.legacy and set(runtime) != set(LEGACY_RUNTIME):
        raise SystemExit(f'ldd found {sorted(runtime)} of the legacy runtime {LEGACY_RUNTIME}')


if __name__ == '__main__':
    main()
