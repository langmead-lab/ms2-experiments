#!/usr/bin/env python3
"""Fetch inputs by logical name into a data directory and verify each one by sha256.

Usage:
  python -m ms2exp.fetch [get] [--data DIR] [--table data/inputs.tsv] (--set toy | --names FILE | NAME ...)
                         [--from-locations L.json ...] [--link] [--mirror DIR_OR_URL] [--verify-only] [--dry-run]
  python -m ms2exp.fetch table [--manifest results/freeze-2885ee8-00e665de/manifest.json] [--out data/inputs.tsv]

get (the default) puts each named input at DIR/<logical name>; DIR is --data, else MS2EXP_DATA.
A file already there with the table's size and sha256 is kept, and a marker under DIR/.verified
records that, so a second run does not hash it again. Any other input is taken, in order, from
the first --from-locations file that maps its name (a {"locations": {name: path}} map such as
data/locations/*.json, copied, or symlinked with --link), from --mirror (a directory or URL prefix
holding files under their table file names), or from the table's URL. Downloads stream to a .part
file, resume with an HTTP Range request after a broken connection, and are renamed into place only
when their sha256 matches; a mismatch removes the .part file and stops with both checksums. Names
with the same sha256 are fetched once and hard-linked (or symlinked with --link). DIR/locations.json
maps every name fetched or verified so far to its path, for ms2exp.catalog hash and link.
--verify-only fetches nothing and fails if any named input is missing or differs.

--set toy names the inputs of the toy run (toy/inputs.txt).

table writes data/inputs.tsv: one row per input of the freeze manifest plus the toy's own inputs,
with its size, sha256, source, deposited file name and URL. Sources are zenodo:<record> (the
three deposits: 23964218 read sets, 23964220 precursors and ropebwt3 FMDs, 23964222 legacy indexes), zenodo:14854401 (Heng Li's human472 record,
unchanged), github:langmead-lab/ms2-experiments@tools-v1 (the release holding the pinned tool
binaries and runtime libraries), rebuild (the large HPRC2 inputs, built from human472.fmd, and the ms2 source bundles), built (made by a build job of
the manifest), recorded (oracle traces a run records) and "deposit to be decided".
"""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import sys
import time
import urllib.error
import urllib.request

from ms2exp.lib import REPO, atomic_json, file_sha256
from ms2exp.paths import data_dir

TABLE = REPO/'data'/'inputs.tsv'
SETS = {'toy': REPO/'toy'/'inputs.txt'}
FIELDS = ('name', 'bytes', 'sha256', 'source', 'file', 'url')
UNDECIDED = 'deposit to be decided'
# The deposits the paper's Data availability note names, by Zenodo record.
DEPOSIT = {'reads': '23964218', 'precursors': '23964220', 'legacy': '23964222'}
DEPOSITED = {'zenodo:' + record for record in DEPOSIT.values()}
ZENODO_URL = 'https://zenodo.org/records/{record}/files/{file}?download=1'
# The pinned ropebwt3, MONI, PHONI and AugPHONI binaries and their runtime libraries are assets
# of one GitHub release of this repository.
RELEASE = 'github:langmead-lab/ms2-experiments@tools-v1'
RELEASE_URL = 'https://github.com/langmead-lab/ms2-experiments/releases/download/tools-v1/{file}'
LI_RECORD = '14854401'
# Inputs of the toy that no freeze manifest has: the mtb152 FMD the mtb152 precursors were made
# from, which the toy's ropebwt3 cells query.
EXTRA = {'ropebwt3/mtb152/index.fmd': dict(
    bytes=19833552, sha256='181c92baeda3f414e50c1e5e87b71ac2c3e54d8bc47237da151aecc192b03381')}


def deposit_file(name):
    """The flat file name an input has in a deposit (Zenodo records have no folders)."""
    return name.replace('/', '__')


def source_of(name, entry):
    top, rest = name.split('/', 1)
    ref = rest.split('/', 1)[0]
    if top == 'reads':
        return 'zenodo:' + DEPOSIT['reads']
    if top in ('moni', 'phoni', 'augphoni'):
        return 'zenodo:' + DEPOSIT['legacy']
    if top == 'precursors':
        return 'rebuild' if ref == 'human472' else 'zenodo:' + DEPOSIT['precursors']
    if top == 'ropebwt3':
        if ref.startswith('chr22-') or name == 'ropebwt3/mtb152/index.fmd':
            return 'zenodo:' + DEPOSIT['precursors']
        if ref == 'human472':
            leaf = rest.split('/', 1)[1]
            return f'zenodo:{LI_RECORD}' if leaf in ('index.fmd.s11.ssa', 'index.fmd.len.gz') else 'rebuild'
        return UNDECIDED
    if top == 'tools':
        return RELEASE
    if top == 'ms2':
        return 'built'
    if top == 'traces':
        return 'recorded'
    if top == 'source':
        return 'rebuild: toy/make_source_bundles.sh from the public ms2 repository'
    raise SystemExit(f'No source rule for {name}')


def built_outputs(results):
    out = {}
    for path in sorted(Path(results).glob('*/outputs.json')):
        out.update(json.loads(path.read_text())['outputs'])
    return out


def cmd_table(args):
    manifest = json.loads(Path(args.manifest).read_text())
    built = built_outputs(Path(args.manifest).parent)
    entries = dict(manifest['inputs'])
    entries.update(EXTRA)
    rows, first = [], {}
    for name in sorted(entries):
        entry = dict(entries[name])
        if entry['sha256'] == 'built' and name in built:
            entry.update(built[name])
        source = source_of(name, entry)
        if name.startswith('tools/ms2-') or name == 'tools/ms_check':
            source = UNDECIDED
        sha = entry['sha256']
        file = url = '-'
        if source in DEPOSITED or source in (RELEASE, UNDECIDED):
            # Byte-identical inputs are deposited once, under the first name that has them.
            file = first.setdefault(sha, deposit_file(name))
            if source == RELEASE:
                url = RELEASE_URL.format(file=file)
            elif source != UNDECIDED:
                url = ZENODO_URL.format(record=source.split(':', 1)[1], file=file)
        elif source == f'zenodo:{LI_RECORD}':
            url = f'https://zenodo.org/records/{LI_RECORD}'
        rows.append(dict(name=name, bytes=entry.get('bytes') or '-', sha256=sha, source=source, file=file, url=url))
    with Path(args.out).open('w', newline='') as f:
        w = csv.DictWriter(f, FIELDS, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    print(f'{args.out}: {len(rows)} inputs')


def load_table(path):
    with Path(path).open(newline='') as f:
        return {row['name']: row for row in csv.DictReader(f, delimiter='\t')}


def wanted(args):
    names = list(args.names_given)
    if args.set:
        names += [line.strip() for line in SETS[args.set].read_text().splitlines()
                  if line.strip() and not line.startswith('#')]
    if args.names:
        names += [line.strip() for line in Path(args.names).read_text().splitlines() if line.strip()]
    if not names:
        raise SystemExit('Name the inputs: --set, --names FILE or NAME ...')
    return list(dict.fromkeys(names))


def verified(data, row, path):
    """True when PATH has the row's size and sha256; a marker remembers a good hash."""
    if not path.exists() or path.stat().st_size != int(row['bytes']):
        return False
    marker = data/'.verified'/row['sha256']
    st = path.stat()
    stamp = f'{path.resolve()} {st.st_size} {st.st_mtime_ns}'
    if marker.exists() and marker.read_text() == stamp:
        return True
    if file_sha256(path) != row['sha256']:
        return False
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(stamp)
    return True


def download(url, dest, row):
    """Stream URL to DEST.part, resuming after a broken connection, then check the sha256."""
    part = dest.with_name(dest.name + '.part')
    size = int(row['bytes'])
    for attempt in range(5):
        have = part.stat().st_size if part.exists() else 0
        if have >= size:
            break
        request = urllib.request.Request(url, headers={'Range': f'bytes={have}-'} if have else {})
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                if have and response.status != 206:
                    have = 0  # the server ignored the range: start again
                with part.open('ab' if have else 'wb') as out:
                    shutil.copyfileobj(response, out, 16 << 20)
        except (urllib.error.URLError, OSError) as e:
            print(f'  {dest.name}: {e}; retrying', file=sys.stderr)
            time.sleep(2 ** attempt)
    actual = file_sha256(part) if part.exists() else 'missing'
    if actual != row['sha256']:
        part.unlink(missing_ok=True)
        raise SystemExit(f"{row['name']}: sha256 {actual}, expected {row['sha256']}")
    part.replace(dest)


def copy_checked(src, dest, row):
    part = dest.with_name(dest.name + '.part')
    shutil.copyfile(src, part)
    actual = file_sha256(part)
    if actual != row['sha256']:
        part.unlink()
        raise SystemExit(f"{row['name']}: {src} has sha256 {actual}, expected {row['sha256']}")
    part.replace(dest)


def public_url(row):
    """True when the row's URL names a file that can be downloaded as is."""
    return row['url'] not in ('', '-') and row['source'] != UNDECIDED and row['source'] != f'zenodo:{LI_RECORD}'


def obtainable(row, args, locations):
    """Where obtain would take the input from, or None when no source has it."""
    if row['name'] in locations:
        return 'location'
    if args.mirror and row['file'] != '-':
        return 'mirror'
    return 'url' if public_url(row) else None


def obtain(row, dest, args, locations):
    """Put the input at DEST from the first source that has it; returns where it came from."""
    if row['name'] in locations:
        src = Path(locations[row['name']])
        if args.link:
            if src.stat().st_size != int(row['bytes']) or file_sha256(src) != row['sha256']:
                raise SystemExit(f"{row['name']}: {src} does not match its sha256")
            dest.symlink_to(src.resolve())
        else:
            copy_checked(src, dest, row)
        return str(src)
    if args.mirror and row['file'] != '-':
        if '://' in args.mirror:
            url = args.mirror.rstrip('/') + '/' + row['file']
            download(url, dest, row)
            return url
        copy_checked(Path(args.mirror)/row['file'], dest, row)
        return str(Path(args.mirror)/row['file'])
    if not public_url(row):
        raise SystemExit(f"{row['name']}: no public copy (source: {row['source']}); "
                         'give --from-locations or --mirror')
    download(row['url'], dest, row)
    return row['url']


def cmd_get(args):
    data = data_dir(args.data)
    table = load_table(args.table)
    names = wanted(args)
    unknown = [n for n in names if n not in table]
    if unknown:
        raise SystemExit('Not in the inputs table: ' + ', '.join(unknown))
    locations = {}
    for path in args.from_locations:
        for name, location in json.loads(Path(path).read_text())['locations'].items():
            locations.setdefault(name, location)
    index_path = data/'locations.json'
    index = json.loads(index_path.read_text())['locations'] if index_path.exists() else {}
    by_sha, missing, unavailable, total = {}, [], [], 0
    for name in names:
        row = table[name]
        if row['bytes'] == '-' or row['sha256'] in ('built', 'pending'):
            raise SystemExit(f'{name} has no fixed content to fetch (sha256 {row["sha256"]})')
        dest = data/name
        if dest.is_symlink() or dest.exists():
            if verified(data, row, dest):
                index[name] = str(dest)
                by_sha.setdefault(row['sha256'], dest)
                print(f'ok      {name}')
                continue
            if args.verify_only:
                missing.append(name + ' (differs)')
                continue
            dest.unlink()
        if args.verify_only:
            missing.append(name)
            continue
        if args.dry_run:
            if obtainable(row, args, locations) is None and row['sha256'] not in by_sha:
                unavailable.append(f"{name} (source: {row['source']})")
                print(f'NONE    {name} (no public copy; source: {row["source"]})')
            else:
                total += int(row['bytes'])
                print(f'fetch   {name} ({int(row["bytes"]):,} bytes)')
            by_sha.setdefault(row['sha256'], dest)
            continue
        total += int(row['bytes'])
        dest.parent.mkdir(parents=True, exist_ok=True)
        twin = by_sha.get(row['sha256'])
        if twin:
            if args.link:
                dest.symlink_to(twin.resolve())
            else:
                os.link(twin, dest)
            print(f'same    {name} (as {twin.relative_to(data)})')
        else:
            origin = obtain(row, dest, args, locations)
            verified(data, row, dest)
            by_sha[row['sha256']] = dest
            print(f'fetched {name} from {origin}')
        index[name] = str(dest)
    if missing:
        raise SystemExit('Missing or different:\n  ' + '\n  '.join(missing))
    if unavailable:
        print(f'{len(names)} inputs; {total:,} bytes to fetch', file=sys.stderr)
        raise SystemExit(f'{len(unavailable)} inputs have no public copy; give --from-locations or --mirror:\n  '
                         + '\n  '.join(unavailable))
    if not args.dry_run:
        atomic_json(index_path, dict(locations=dict(sorted(index.items()))))
    print(f'{len(names)} inputs in {data}' + (f'; {total:,} bytes to fetch' if args.dry_run else ''))


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == 'table':
        p = argparse.ArgumentParser(prog='ms2exp.fetch table')
        p.add_argument('--manifest', default=str(REPO/'results/freeze-2885ee8-00e665de/manifest.json'))
        p.add_argument('--out', default=str(TABLE))
        return cmd_table(p.parse_args(argv[1:]))
    if argv and argv[0] == 'get':
        argv = argv[1:]
    p = argparse.ArgumentParser(prog='ms2exp.fetch', description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('names_given', nargs='*', metavar='NAME')
    p.add_argument('--data')
    p.add_argument('--table', default=str(TABLE))
    p.add_argument('--set', choices=sorted(SETS))
    p.add_argument('--names', help='file of logical names, one per line')
    p.add_argument('--from-locations', nargs='+', default=[], metavar='L.json')
    p.add_argument('--link', action='store_true', help='symlink --from-locations files instead of copying')
    p.add_argument('--mirror', help='directory or URL prefix holding files under their table file names')
    p.add_argument('--verify-only', action='store_true')
    p.add_argument('--dry-run', action='store_true')
    return cmd_get(p.parse_args(argv))


if __name__ == '__main__':
    main()
