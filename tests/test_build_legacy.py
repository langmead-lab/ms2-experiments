"""The legacy-tool build recipe, ms2exp.build.legacy."""
import json

from ms2exp.build import legacy
from ms2exp.lib import REPO

LOCK = json.loads((REPO/'tools.lock.json').read_text())


def dry_run(tmp_path, *extra):
    args = legacy.parse(['--dry-run', '--work', str(tmp_path/'w'), '--out', str(tmp_path/'o'), *extra])
    return [' '.join(c['argv']) for c in legacy.Builder(args, LOCK).main()['commands']]


def test_lock_pins_every_fetched_dependency_by_commit():
    fetch = LOCK['legacy_build']['fetchcontent']
    assert set(fetch) == {'moni', 'phoni', 'augphoni'}
    for tool, deps in fetch.items():
        assert {'SDSL', 'R-INDEX', 'SHAPED_SLP', 'BWT2LCP'} <= set(deps), tool
        for dep in deps.values():
            assert len(dep['commit']) == 40 and dep['repo'].startswith('https://')
    for patches in LOCK['legacy_build']['dependency_patches'].values():
        for patch in patches:
            assert (REPO/patch['patch']).is_file()
    for tool in ('moni', 'phoni', 'augphoni'):
        for patch in LOCK['legacy_build']['tools'][tool]['patches']:
            assert (REPO/patch).is_file()


def test_cmake_tools_build_offline_from_pinned_dependencies(tmp_path):
    commands = dry_run(tmp_path, '--tools', 'phoni')
    configure = next(c for c in commands if ' -S ' in c)
    assert '-DFETCHCONTENT_FULLY_DISCONNECTED=ON' in configure
    assert f"-DFETCHCONTENT_SOURCE_DIR_SDSL={tmp_path/'w'/'phoni-deps'/'sdsl-src'}" in configure
    assert '-ftrivial-auto-var-init=zero' in configure and 'prefix-map' not in configure
    assert any(' am ' in c and c.endswith('patches/phoni.patch') for c in commands)
    assert any(c.endswith('patches/sux-utility.patch') for c in commands)
    assert any(c.endswith('patches/bwt2lcp-cstdint.patch') for c in commands)
    build = next(c for c in commands if ' --build ' in c)
    assert ' phoni ' in build and build.endswith('-j2')


def test_freeze_paths_map_source_build_and_dependency_paths(tmp_path):
    commands = dry_run(tmp_path, '--tools', 'moni,ropebwt3', '--freeze-paths', '--jobs', '1')
    configure = next(c for c in commands if ' -S ' in c)
    w = tmp_path/'w'
    for new, old in ((w/'moni-src', '/scratch16/blangme2/ms-bakery/moni'),
                     (w/'moni-build', '/scratch16/blangme2/ms-bakery/moni/build-bakery'),
                     (w/'moni-deps', '/scratch16/blangme2/ms-bakery/moni/build-bakery/_deps')):
        assert f'-fmacro-prefix-map={new}={old}' in configure
    make = next(c for c in commands if c.startswith('make '))
    # The Makefile's own CFLAGS and CPPFLAGS (with OpenMP) stay in force.
    assert 'CFLAGS' not in make and 'CPPFLAGS' not in make
    assert f"INCLUDES=-ffile-prefix-map={w/'ropebwt3'}=/scratch16/blangme2/ms-bakery/ropebwt3" in make
    assert make.startswith('make -j1 ')
