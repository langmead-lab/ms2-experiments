#!/usr/bin/env python3
"""Download Guccione et al. (PRJEB83637) skin metagenome runs from ENA and check each file's size and MD5.

The run table is manifests/guccione-ena-runs.tsv (the ENA read_run report of the study, 127 runs). By
default only ERR14066769, the run of the skin read set, is fetched; --all fetches every skin run, which
ms2exp.prep.skin_screen needs. Files land in ROOT/skin/RUN/ and ROOT/skin-verified.json lists them with
byte count, MD5 and sha256. A file already present is checked and reused.

Usage: python3 data/download/download_skin.py ROOT [--all | --runs RUN ...]
"""
import argparse
import concurrent.futures
import csv
import hashlib
import json
from pathlib import Path
import subprocess

TABLE = Path(__file__).resolve().parent/'manifests/guccione-ena-runs.tsv'
SET_RUN = 'ERR14066769'


def files(table=TABLE, runs=None):
    rows = list(csv.DictReader(Path(table).open(), delimiter='\t'))
    out = []
    for r in rows:
        if r['scientific_name'] != 'human skin metagenome' or (runs is not None and r['run_accession'] not in runs):
            continue
        for url, md5, n in zip(r['fastq_ftp'].split(';'), r['fastq_md5'].split(';'), r['fastq_bytes'].split(';')):
            out.append(dict(cohort='skin', run=r['run_accession'], url='https://' + url, md5=md5, bytes=int(n),
                            relative_path='skin/' + r['run_accession'] + '/' + url.rsplit('/', 1)[-1]))
    if runs is not None and {f['run'] for f in out} != set(runs):
        raise ValueError('Run not in the skin table: ' + ', '.join(sorted(set(runs) - {f['run'] for f in out})))
    return out


def digest(path):
    m, s = hashlib.md5(), hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(4 << 20), b''):
            m.update(block)
            s.update(block)
    return m.hexdigest(), s.hexdigest()


def get(root, f):
    dest = root/f['relative_path']
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        m, s = digest(dest)
        if dest.stat().st_size == f['bytes'] and m == f['md5']:
            return dict(f, sha256=s, reused=True)
        raise RuntimeError('Existing file does not match ENA: ' + str(dest))
    part = dest.with_suffix(dest.suffix + '.part')
    subprocess.run(['curl', '--fail', '--location', '--silent', '--show-error', '--retry', '4', '--retry-delay', '5',
                    '--connect-timeout', '30', '--max-time', '7200', '--continue-at', '-', '--output', str(part),
                    f['url']], check=True)
    if part.stat().st_size != f['bytes']:
        raise RuntimeError('Size mismatch: ' + str(part))
    m, s = digest(part)
    if m != f['md5']:
        raise RuntimeError('MD5 mismatch: ' + str(part))
    part.replace(dest)
    return dict(f, sha256=s, reused=False)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('root', type=Path)
    g = p.add_mutually_exclusive_group()
    g.add_argument('--all', action='store_true')
    g.add_argument('--runs', nargs='+')
    a = p.parse_args(argv)
    todo = files(runs=None if a.all else (a.runs or [SET_RUN]))
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        verified = list(pool.map(lambda f: get(a.root, f), todo))
    (a.root/'skin-verified.json').write_text(json.dumps(verified, indent=2) + '\n')
    print(len(verified), 'files verified', flush=True)


if __name__ == '__main__':
    main()
