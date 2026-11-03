#!/usr/bin/env python3
"""Build ropebwt3, MONI, PHONI and Augmented PHONI from their public commits, as the freeze did.

Usage:
  ms2exp/build/build_legacy.sh [--tools ropebwt3,moni,phoni,augphoni] [--lock tools.lock.json]
      [--work DIR] [--out DIR] [--cc gcc] [--cxx g++] [--cmake cmake] [--jobs N]
      [--freeze-paths] [--static-runtime] [--dry-run]

Everything comes from tools.lock.json: each tool's repository and commit, the patches in patches/,
and the legacy_build block, which holds the CMake flags and targets of the freeze builds and the
commit of every dependency that CMake's FetchContent fetched at configure time. Several of those
dependencies are named upstream by branch (master, develop, main), so this script clones each one
at its recorded commit itself and points CMake at the clone with FETCHCONTENT_SOURCE_DIR_<NAME>,
with FETCHCONTENT_FULLY_DISCONNECTED so that CMake fetches nothing. Two dependencies get a one-line
include that GCC 13 needs (legacy_build.dependency_patches): bwt2lcp and sux, a submodule of ShapedSlp.

PHONI and Augmented PHONI are built from their upstream base commits with patches/phoni.patch and
patches/aug_phoni.patch, which hold the commits of the upstream pull requests. The freeze binaries
were built from local branches whose code is the same and whose headers differ only in comments.

The freeze builds used GCC 13.1.0 and CMake 3.27.7 (legacy_build.compiler and legacy_build.cmake).
Another GCC 13 should give the same machine code; GCC 12 or later is required, since PHONI and
Augmented PHONI need -ftrivial-auto-var-init=zero.

--freeze-paths adds prefix maps that rewrite the source, build and dependency paths of this build to
the freeze build's paths (legacy_build.tools.*.freeze_paths): -fmacro-prefix-map for the CMake tools,
whose binaries embed source paths through __FILE__, and -ffile-prefix-map for ropebwt3, whose
binary also carries debug information. With the freeze compiler this aims at byte-identical
binaries; tools/build-legacy.json says whether each one matched its pinned sha256.

--static-runtime links libstdc++ and libgcc of MONI, PHONI and Augmented PHONI statically, for hosts
without a GCC 13 runtime. Those binaries differ from the pinned ones by construction. ropebwt3 is C
and needs only the OpenMP runtime (libgomp), which this option leaves dynamic.

Writes OUT/tools/<tool>/<binary> and OUT/tools/build-legacy.json (commands, compiler and CMake
versions, dependency commits, flags and the sha256 of every binary).
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from ms2exp.lib import REPO, atomic_json, file_sha256

TOOLS = ('ropebwt3', 'moni', 'phoni', 'augphoni')
OUT_NAME = {'ropebwt3': 'ropebwt3', 'moni': 'moni-ms', 'phoni': 'phoni', 'augphoni': 'aug_phoni'}
GIT_ID = ['-c', 'user.name=build_legacy', '-c', 'user.email=build_legacy@localhost']


class Builder:
    def __init__(self, args, lock):
        self.args, self.lock = args, lock
        self.legacy = lock['legacy_build']
        self.work = Path(args.work).resolve()
        self.out = Path(args.out).resolve()
        self.record = dict(started=time.time(), argv=sys.argv, lock=str(Path(args.lock).resolve()),
                           freeze_paths=args.freeze_paths, static_runtime=args.static_runtime, tools={},
                           commands=[])

    def run(self, argv, cwd=None, env=None, check=True):
        argv = [str(a) for a in argv]
        print('+', ' '.join(argv), flush=True)
        if self.args.dry_run:
            self.record['commands'].append(dict(argv=argv, cwd=str(cwd) if cwd else None))
            return ''
        started = time.time()
        result = subprocess.run(argv, cwd=cwd, env=env, capture_output=True, text=True)
        self.record['commands'].append(dict(argv=argv, cwd=str(cwd) if cwd else None, returncode=result.returncode,
                                            seconds=round(time.time() - started, 3)))
        if check and result.returncode:
            sys.stdout.write(result.stdout[-4000:] + result.stderr[-4000:])
            raise SystemExit(f'Command failed ({result.returncode}): {" ".join(argv)}')
        return result.stdout

    def clone(self, repo, commit, dest, submodules=None):
        """dest at commit; submodules maps path to the commit the dependency records, checked after update."""
        if (dest/'.git').exists() and self.run(['git', 'rev-parse', 'HEAD'], dest).strip() == commit:
            return
        if dest.exists() and not self.args.dry_run:
            raise SystemExit(f'{dest} exists at another commit; give a fresh --work')
        self.run(['git', 'clone', '-q', repo, dest])
        self.run(['git', '-c', 'advice.detachedHead=false', 'checkout', '-q', commit], dest)
        if submodules is not None:
            self.run(['git', 'submodule', 'update', '-q', '--init', '--recursive'], dest)
            for path, want in submodules.items():
                head = self.run(['git', 'rev-parse', 'HEAD'], dest/path).strip()
                if head != want and not self.args.dry_run:
                    raise SystemExit(f'{dest/path} is at {head}, not {want}')
        head = self.run(['git', 'rev-parse', 'HEAD'], dest).strip()
        if head != commit and not self.args.dry_run:
            raise SystemExit(f'{dest} is at {head}, not {commit}')

    def compiler(self):
        version = self.run([self.args.cxx, '-dumpfullversion']).strip()
        cmake = self.run([self.args.cmake, '--version']).splitlines()
        self.record['compiler'] = dict(cc=self.args.cc, cxx=self.args.cxx, version=version,
                                       cc_version=self.run([self.args.cc, '--version']).splitlines()[:1],
                                       cxx_version=self.run([self.args.cxx, '--version']).splitlines()[:1],
                                       cmake=cmake[:1])
        if self.args.dry_run:
            return
        want = self.legacy['compiler']['version']
        if version != want:
            print(f'Warning: {self.args.cxx} is GCC {version}; the freeze binaries were built with {want}', flush=True)
        if int(version.split('.')[0]) < 12 and {'phoni', 'augphoni'} & set(self.args.tools):
            raise SystemExit('PHONI and Augmented PHONI need GCC 12 or later for -ftrivial-auto-var-init=zero')

    def install(self, tool, built, info):
        dest = self.out/'tools'/tool/OUT_NAME[tool]
        info['binary'] = str(dest)
        if self.args.dry_run:
            return
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(built, dest)
        dest.chmod(0o755)
        info['sha256'] = file_sha256(dest)
        info['pinned_sha256'] = self.lock[tool]['binary_sha256']
        info['matches_freeze'] = info['sha256'] == info['pinned_sha256']
        print(f"{tool}: {info['sha256']} ({'matches' if info['matches_freeze'] else 'differs from'} the pin)",
              flush=True)

    def ropebwt3(self):
        pin, spec = self.lock['ropebwt3'], self.legacy['tools']['ropebwt3']
        src = self.work/'ropebwt3'
        self.clone(pin['repo'], pin['commit'], src)
        # Command-line CFLAGS or CPPFLAGS would replace the Makefile's OpenMP additions, so the prefix map
        # goes through INCLUDES, which the Makefile leaves empty.
        argv = ['make', f"-j{self.args.jobs or spec['jobs']}", f'CC={self.args.cc}', *spec['make_args']]
        includes = ''
        if self.args.freeze_paths:
            includes = f"-ffile-prefix-map={src}={spec['freeze_paths']['source']}"
            argv.append(f'INCLUDES={includes}')
        self.run(argv, src)
        info = dict(repo=pin['repo'], commit=pin['commit'], make=argv, includes=includes, patches=[])
        self.install('ropebwt3', src/spec['binary'], info)
        self.record['tools']['ropebwt3'] = info

    def cmake_tool(self, tool):
        pin, spec = self.lock[tool], self.legacy['tools'][tool]
        src, build, deps = (self.work/f'{tool}-src', self.work/f'{tool}-build', self.work/f'{tool}-deps')
        self.clone(pin['repo'], pin['commit'], src)
        patches = [REPO/p for p in spec.get('patches', [])]
        for patch in patches:
            if patch.read_text().startswith('From '):
                self.run(['git', *GIT_ID, 'am', '-q', '--committer-date-is-author-date', patch], src)
            else:
                self.run(['git', 'apply', patch], src)
        info = dict(repo=pin['repo'], commit=pin['commit'], patches=[], dependencies={})
        for patch in patches:
            info['patches'].append(dict(path=str(patch.relative_to(REPO)), sha256=file_sha256(patch)))
        defines = []
        for name, dep in sorted(self.legacy['fetchcontent'][tool].items()):
            dest = deps/f'{name.lower()}-src'
            self.clone(dep['repo'], dep['commit'], dest, submodules=dep.get('submodules', {}))
            patches = self.legacy['dependency_patches'].get(name, [])
            for patch in patches:
                self.run(['git', 'apply', REPO/patch['patch']], dest/patch['dir'])
            info['dependencies'][name] = dict(dep, patches=patches)
            defines.append(f'-DFETCHCONTENT_SOURCE_DIR_{name}={dest}')
        cxx_flags, c_flags = spec['cxx_flags'], ''
        if self.args.freeze_paths:
            freeze = spec['freeze_paths']
            maps = ' '.join(f'-fmacro-prefix-map={new}={old}' for new, old in
                            ((src, freeze['source']), (build, freeze['build']), (deps, freeze['deps'])))
            cxx_flags, c_flags = f'{cxx_flags} {maps}', maps
        argv = [self.args.cmake, '-S', src, '-B', build, f"-DCMAKE_BUILD_TYPE={spec['build_type']}",
                f'-DCMAKE_C_COMPILER={self.args.cc}', f'-DCMAKE_CXX_COMPILER={self.args.cxx}',
                f'-DCMAKE_CXX_FLAGS={cxx_flags}', f'-DCMAKE_C_FLAGS={c_flags}',
                '-DFETCHCONTENT_FULLY_DISCONNECTED=ON', *defines]
        if self.args.static_runtime:
            argv.append('-DCMAKE_EXE_LINKER_FLAGS=-static-libstdc++ -static-libgcc')
        self.run(argv)
        self.run([self.args.cmake, '--build', build, '--target', *spec['targets'],
                  f"-j{self.args.jobs or spec['jobs']}"])
        info.update(cxx_flags=cxx_flags, c_flags=c_flags, targets=spec['targets'])
        self.install(tool, build/spec['binary'], info)
        self.record['tools'][tool] = info

    def main(self):
        if not self.args.dry_run:
            self.work.mkdir(parents=True, exist_ok=True)
        self.compiler()
        for tool in self.args.tools:
            self.ropebwt3() if tool == 'ropebwt3' else self.cmake_tool(tool)
        self.record['finished'] = time.time()
        if not self.args.dry_run:
            atomic_json(self.out/'tools'/'build-legacy.json', self.record)
        return self.record


def parse(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--tools', default=','.join(TOOLS), help='comma-separated subset of ' + ', '.join(TOOLS))
    p.add_argument('--lock', default=str(REPO/'tools.lock.json'))
    p.add_argument('--work', default='legacy-build', help='checkouts and build trees (must be fresh)')
    p.add_argument('--out', default='legacy-out', help='writes OUT/tools/<tool>/<binary> and build-legacy.json')
    p.add_argument('--cc', default=os.environ.get('CC', 'gcc'))
    p.add_argument('--cxx', default=os.environ.get('CXX', 'g++'))
    p.add_argument('--cmake', default='cmake')
    p.add_argument('--jobs', type=int, help='parallel build jobs (default: the freeze builds\' 2)')
    p.add_argument('--freeze-paths', action='store_true', help='map build paths onto the freeze paths')
    p.add_argument('--static-runtime', action='store_true', help='link libstdc++ and libgcc statically')
    p.add_argument('--dry-run', action='store_true', help='print the commands only')
    args = p.parse_args(argv)
    args.tools = [t for t in args.tools.split(',') if t]
    unknown = set(args.tools) - set(TOOLS)
    if unknown:
        p.error(f'unknown tools: {sorted(unknown)}')
    return args


def main(argv=None):
    args = parse(argv)
    lock = json.loads(Path(args.lock).read_text())
    Builder(args, lock).main()


if __name__ == '__main__':
    main()
