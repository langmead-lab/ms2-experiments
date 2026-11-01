#!/usr/bin/env python3
"""Print (or, with --submit, run) the sbatch commands of the freeze jobs on Rockfish.

Usage:
  python -m launchers.rockfish.submit --manifest M --inputs-root DIR --results-root DIR [--jobs JOB ...] [--label orbitms]
            [--time HH:MM:SS] [--submit]

Each job goes to a partition sized from its largest index: bigmem (account blangme2_bigmem,
qos_bigmem, 1024 GB) when any index of the job is over 100 GiB, else parallel (46 CPUs, 184 GB).
Every job is exclusive, so its timings see one job per node, as on AWS (D5). Jobs whose done
marker exists in the results tree are skipped. Without --submit nothing is submitted; check the
printed commands and the resource caps (rockfish skill: slurm_preflight.sh) first.
"""
import argparse
from pathlib import Path
import shlex
import subprocess

from ms2exp.lib import load_manifest

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

BIG = 100 << 30
PARTITIONS = {
    'bigmem': ['--exclusive', '--partition=bigmem', '--account=blangme2_bigmem', '--qos=qos_bigmem', '--mem=1000G',
               '--cpus-per-task=48'],
    'parallel': ['--exclusive', '--partition=parallel', '--account=blangme2', '--mem=184G', '--cpus-per-task=46'],
}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest', required=True)
    p.add_argument('--inputs-root', required=True)
    p.add_argument('--results-root', required=True)
    p.add_argument('--jobs', nargs='*')
    p.add_argument('--label', default='orbitms')
    p.add_argument('--time', default='24:00:00')
    p.add_argument('--partition', choices=sorted(PARTITIONS), help='override the size rule')
    p.add_argument('--logs', default=None, help='log directory (default RESULTS_ROOT/MANIFEST_ID/logs)')
    p.add_argument('--submit', action='store_true')
    args = p.parse_args()
    manifest = load_manifest(args.manifest)
    results = Path(args.results_root)
    logs = Path(args.logs or results/manifest['id']/'logs')
    for job in manifest['jobs']:
        if job.get('kind') == 'build' or 'rockfish' not in job['platforms'] or (args.jobs and job['id'] not in args.jobs):
            continue
        if (results/manifest['id']/job['id']/'job.done').exists():
            print('# done:', job['id'])
            continue
        largest = max(((manifest['inputs'][n].get('bytes') or 0) for n in job['inputs']
                       if n.startswith(('ms2/', 'ropebwt3/'))), default=0)
        partition = args.partition or ('bigmem' if largest > BIG else 'parallel')
        argv = ['sbatch', '-J', f"{args.label}-freeze-{job['id']}", f'--time={args.time}',
                '-o', str(logs/f"{job['id']}-%j.log"), *PARTITIONS[partition], f'--export=ALL,FREEZE_HOME={REPO}',
                str(HERE/'run_job.sbatch'), str(Path(args.manifest).resolve()), job['id'],
                args.inputs_root, args.results_root]
        print(shlex.join(argv))
        if args.submit:
            logs.mkdir(parents=True, exist_ok=True)
            subprocess.run(argv, check=True)


if __name__ == '__main__':
    main()
