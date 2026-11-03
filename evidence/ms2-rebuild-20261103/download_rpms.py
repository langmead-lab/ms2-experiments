# Resolves and downloads packages with the dnf API, which needs no root for downloads.
import shutil, sys
import dnf
base = dnf.Base()
base.conf.cachedir = '/tmp/dnfc'
base.conf.substitutions.update_from_etc('/')
base.read_all_repos()
base.fill_sack(load_system_repo=True)
for name in sys.argv[1:]:
    base.install(name)
base.resolve()
pkgs = list(base.transaction.install_set)
base.download_packages(pkgs)
for p in pkgs:
    shutil.copy(p.localPkg(), '/rpms/')
    print(p.name, p.evr, p.reponame)
