#!/usr/bin/env python3
"""Derive a rerun manifest from an earlier run of the same matrix and seed its results.

Usage: seed_rerun.py OLD_MANIFEST NEW_MANIFEST --results-root ROOT [--note TEXT] [--dry-run]

The rerun keeps the old manifest's jobs, cells and inputs and takes the protocol from
make_manifest.PROTOCOL, so a protocol change (for example transparent_hugepages) gives it a new id
and its own results prefix. Into that prefix it copies from the old run's prefix:

- each build job's job.done and outputs.json. Built files live under inputs/ by logical name, so
  the rerun's jobs stage the same bytes and no build runs again;
- the recorded oracle traces (traces/);
- every untimed job, whole. Its cells (operation counts, A1) do not depend on the timing protocol.

Timed jobs start empty, so every timed cell runs again. It refuses when an old build job has no
job.done. The new prefix records what it copied in seeded-from.json. ROOT is a directory or an S3
prefix (config.json's bench_s3_root plus /results); S3 copies run server side.
"""
import argparse
import copy
import json
from pathlib import Path

from ms2exp.lib import Store, json_sha256, load_manifest, now
from ms2exp.manifest import PROTOCOL


def derive(old, note):
    new = copy.deepcopy(old)
    del new['id']
    new.update(created=now(), protocol=PROTOCOL, rerun_of=old['id'], note=note or old.get('note'))
    new['id'] = f"{new['matrix']}-{new['ms2']['commit'][:7]}-{json_sha256(new)[:8]}"
    return new


def copy_names(store, names, src, dest, dry):
    for name in names:
        target = dest + name[len(src):]
        if dry:
            continue
        if store.s3:
            store._aws('s3', 'cp', '--only-show-errors', store.url(name), store.url(target))
        else:
            store.put(store.local_path(name), target)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('old_manifest')
    p.add_argument('new_manifest')
    p.add_argument('--results-root', required=True)
    p.add_argument('--note')
    p.add_argument('--dry-run', action='store_true', help='write the manifest and list the copies only')
    args = p.parse_args()
    old = load_manifest(args.old_manifest)
    new = derive(old, args.note)
    store = Store(args.results_root)
    present = set(store.names(old['id']))
    jobs = [j for j in old['jobs'] if 'aws' in j['platforms']]
    builds = [j['id'] for j in jobs if j['kind'] == 'build']
    missing = [b for b in builds if f"{old['id']}/{b}/job.done" not in present]
    if missing:
        raise SystemExit('Build jobs without job.done in the old run: ' + ', '.join(missing))
    picked = {}
    for b in builds:
        picked[b] = [f"{old['id']}/{b}/{n}" for n in ('job.done', 'outputs.json')]
    for j in jobs:
        if j.get('pool') == 'untimed':
            picked[j['id']] = sorted(n for n in present if n.startswith(f"{old['id']}/{j['id']}/"))
    picked['traces'] = sorted(n for n in present if n.startswith(f"{old['id']}/traces/"))
    Path(args.new_manifest).write_text(json.dumps(new, indent=1) + '\n')
    total = 0
    for what, names in picked.items():
        copy_names(store, names, old['id'], new['id'], args.dry_run)
        total += len(names)
        print(f'{what}: {len(names)} files')
    record = dict(rerun_of=old['id'], seeded=now(), protocol_change={
        k: (old['protocol'].get(k), new['protocol'].get(k))
        for k in set(old['protocol']) | set(new['protocol']) if old['protocol'].get(k) != new['protocol'].get(k)},
        copied={what: len(names) for what, names in picked.items()})
    if not args.dry_run:
        store.put_json(record, f"{new['id']}/seeded-from.json")
    timed = sum(1 for j in jobs if j.get('pool') == 'timed')
    print(f"{new['id']}: {total} files {'to copy' if args.dry_run else 'copied'} from {old['id']}; "
          f'{timed} timed jobs start empty')


if __name__ == '__main__':
    main()
