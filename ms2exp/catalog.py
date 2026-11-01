#!/usr/bin/env python3
"""Inputs catalog of the freeze benchmark: logical names, files, checksums, and the inputs tree.

A locations file maps each logical name to a file on Rockfish:
  {"locations": {"reads/human472/human_short-ms-60mb-acgt-v1/reads.fa": "/scratch16/...", ...}}
Entries whose value is a directory contribute every file below it, under the logical name
plus the relative path (for tools/ trees written by build_tools.py).

Commands:
  catalog.py hash --locations L.json [L2.json ...] --out inputs.json
      sha256 and size of every location. Unchanged files (same path, size and mtime as in an
      existing --out) keep their recorded sha256, so rehashing after a new build is cheap.
  catalog.py link --locations L.json [...] --root DIR
      Make DIR/<logical name> a symlink to each location: the --inputs-root of Rockfish runs.
  catalog.py publish --manifest M --root DIR --s3 s3://bucket/prefix [--dry-run]
      Upload every input of the manifest from DIR to the S3 prefix under its logical name,
      checking the sha256 first; objects already present with the same size are skipped.
      Traces recorded by Rockfish runs (results traces/) are published with --traces DIR. Inputs an
      all-AWS manifest builds (sha256 "built") are skipped: their build jobs publish them.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from ms2exp.lib import Store, atomic_json, file_sha256, load_manifest


def locations(paths):
    out = {}
    for path in paths:
        for name, location in json.loads(Path(path).read_text())['locations'].items():
            location = Path(location)
            if location.is_dir():
                for f in sorted(p for p in location.rglob('*') if p.is_file()):
                    out[f'{name}/{f.relative_to(location)}'] = f
            else:
                out[name] = location
    return out


def cmd_hash(args):
    previous = {}
    if Path(args.out).exists():
        previous = json.loads(Path(args.out).read_text())['inputs']
    inputs = {}
    for name, path in sorted(locations(args.locations).items()):
        st = path.stat()
        old = previous.get(name, {})
        if (old.get('path'), old.get('bytes'), old.get('mtime_ns')) == (str(path), st.st_size, st.st_mtime_ns):
            inputs[name] = old
            continue
        print('hash', name, file=sys.stderr, flush=True)
        inputs[name] = dict(sha256=file_sha256(path), bytes=st.st_size, path=str(path), mtime_ns=st.st_mtime_ns)
    atomic_json(args.out, dict(format='ms2-freeze-inputs-v1', inputs=inputs))
    print(f'{len(inputs)} inputs in {args.out}')


def cmd_link(args):
    root = Path(args.root)
    for name, path in sorted(locations(args.locations).items()):
        link = root/name
        link.parent.mkdir(parents=True, exist_ok=True)
        if link.is_symlink() or link.exists():
            if link.is_symlink() and Path(os.readlink(link)) == path.resolve():
                continue
            link.unlink()
        os.symlink(path.resolve(), link)
    print(f'linked under {root}')


def cmd_publish(args):
    manifest = load_manifest(args.manifest)
    root, s3 = Path(args.root), Store(args.s3)
    for name, entry in manifest['inputs'].items():
        if entry['sha256'] == 'built':
            continue  # an AWS build job publishes it (build_job.py)
        if entry['sha256'] == 'pending':
            source = Path(args.traces)/name if args.traces and name.startswith('traces/') else None
            if source is None or not source.exists():
                print('pending, not published:', name)
                continue
        else:
            source = root/name
        listing = subprocess.run(['aws', 's3', 'ls', s3.url(name)], capture_output=True, text=True).stdout.split()
        if len(listing) >= 3 and int(listing[2]) == source.stat().st_size:
            continue
        if entry['sha256'] != 'pending' and file_sha256(source) != entry['sha256']:
            raise SystemExit('Checksum differs for ' + name)
        print('upload', name, flush=True)
        if not args.dry_run:
            s3.put(source, name)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    h = sub.add_parser('hash')
    h.add_argument('--locations', nargs='+', required=True)
    h.add_argument('--out', required=True)
    l = sub.add_parser('link')
    l.add_argument('--locations', nargs='+', required=True)
    l.add_argument('--root', required=True)
    u = sub.add_parser('publish')
    u.add_argument('--manifest', required=True)
    u.add_argument('--root', required=True)
    u.add_argument('--s3', required=True)
    u.add_argument('--traces')
    u.add_argument('--dry-run', action='store_true')
    args = p.parse_args()
    dict(hash=cmd_hash, link=cmd_link, publish=cmd_publish)[args.command](args)


if __name__ == '__main__':
    main()
