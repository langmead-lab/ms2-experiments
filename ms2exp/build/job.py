#!/usr/bin/env python3
"""Run one build job of an all-AWS freeze manifest and publish what it builds as inputs.

Usage:
  python -m ms2exp.build.job JOB_ID --manifest MANIFEST --inputs-root ROOT --results-root ROOT --stage DIR [--evict]

A build job (ms2exp.manifest --all-aws) lists pinned inputs, outputs by logical name, and steps:
argv lists whose tokens name staged inputs ({in:NAME}), outputs ({out:NAME}), the output
directory ({outroot}), a scratch directory ({work}), the CPU count ({threads}) and the freeze
code directory ({code}). The steps run in order under /usr/bin/time, without NUMA binding (a
build is not timed for the paper and may need both sockets' memory). A step whose done marker
exists in the work directory is skipped, so a rerun on the same box resumes after the last
finished step.

When every step has run, each output is hashed and uploaded to the inputs store under its logical
name, and RESULTS_ROOT/MANIFEST_ID/JOB_ID/outputs.json records its sha256 and size. Jobs that use
an output carry sha256 "built" in the manifest; the runner takes the checksum from outputs.json
and verifies the staged copy against it, and worker.py offers those jobs only once job.done of
the build exists. A rerun of a finished build job does nothing.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import time

from ms2exp.lib import REPO, Store, atomic_json, file_sha256
from ms2exp.run.runner import Runner, common, log

STEP_TIMEOUT = 48 * 3600


class Builder(Runner):
    def expand(self, token, out_root, work):
        if token.startswith('{in:') and token.endswith('}'):
            return str(self.staged[token[4:-1]])
        if token.startswith('{out:') and token.endswith('}'):
            path = out_root/token[5:-1]
            path.parent.mkdir(parents=True, exist_ok=True)
            return str(path)
        return (token.replace('{outroot}', str(out_root)).replace('{work}', str(work))
                .replace('{threads}', str(os.cpu_count())).replace('{code}', str(REPO)))

    def build(self):
        job = self.job
        if job.get('kind') != 'build':
            raise SystemExit(f"{job['id']} is not a build job")
        record_name = f'{self.prefix}/outputs.json'
        if self.results.exists(f'{self.prefix}/job.done'):
            log(job['id'], 'already built')
            return 0
        attempt = dict(job=job['id'], manifest_id=self.manifest['id'], host=self.host, started=time.time(),
                       argv=sys.argv, steps=[])
        attempt['staging'] = self.stage_names(set(job['inputs']))
        out_root = self.work/'outputs'
        scratch = self.work/'scratch'
        scratch.mkdir(parents=True, exist_ok=True)
        status = 0
        for number, step in enumerate(job['steps'], 1):
            label = f"{number:02d}-{step['label']}"
            marker = self.work/'steps'/(label + '.done')
            if marker.exists():
                attempt['steps'].append(json.loads(marker.read_text()))
                log(job['id'], label, 'done earlier')
                continue
            argv = [self.expand(t, out_root, scratch) for t in step['argv']]
            log(job['id'], label, ' '.join(argv))
            folder = self.work/'steps'/label
            folder.mkdir(parents=True, exist_ok=True)
            try:
                row = common.execute(argv, folder, label, timeout=step.get('timeout_s', STEP_TIMEOUT))
            except RuntimeError as error:
                attempt['steps'].append(dict(label=label, argv=argv, error=str(error)))
                status = 1
                break
            done = dict(label=label, argv=row['argv'], wall_seconds=row['wall_seconds'],
                        peak_rss_kib=row.get('peak_rss_kib'), unverified=step.get('unverified', False))
            atomic_json(marker, done)
            attempt['steps'].append(done)
            self.results.put(folder/(label + '.stderr'), f'{self.prefix}/steps/{label}.stderr')
        if status == 0:
            outputs = {}
            inputs = Store(self.args.inputs_root)
            for name in job['outputs']:
                path = out_root/name
                if not path.is_file():
                    attempt['missing_output'] = name
                    status = 1
                    break
                outputs[name] = dict(sha256=file_sha256(path), bytes=path.stat().st_size)
                log(job['id'], 'publish', name, outputs[name]['sha256'][:12])
                inputs.put(path, name)
            if status == 0:
                self.results.put_json(dict(job=job['id'], manifest_id=self.manifest['id'], outputs=outputs,
                                           host=self.host, finished=time.time()), record_name)
                self.results.put_json(dict(job=job['id'], finished=time.time(), outputs=len(outputs)),
                                      f'{self.prefix}/job.done')
                # Later jobs on this box stage the outputs from the cache instead of from S3.
                for name, out in outputs.items():
                    folder = self.cache/out['sha256']
                    folder.mkdir(parents=True, exist_ok=True)
                    shutil.move(str(out_root/name), folder/Path(name).name)
                    (folder/'verified').write_text(json.dumps(dict(name=name, sha256=out['sha256'],
                                                                   at=time.time())) + '\n')
                shutil.rmtree(self.work, ignore_errors=True)
        attempt.update(finished=time.time(), status=status)
        stamp = time.strftime('%Y%m%dT%H%M%S')
        self.results.put_json(attempt, f"{self.prefix}/attempts/{stamp}-{self.host['hostname']}.json")
        return status


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('job_id')
    p.add_argument('--manifest', required=True)
    p.add_argument('--inputs-root', required=True)
    p.add_argument('--results-root', required=True)
    p.add_argument('--stage', required=True)
    p.add_argument('--platform', default='aws', choices=('aws', 'rockfish'))
    p.add_argument('--evict', action='store_true')
    args = p.parse_args()
    # Runner options that only cell jobs use.
    args.no_numa, args.record_traces, args.keep_work = True, False, False
    sys.exit(Builder(args).build())


if __name__ == '__main__':
    main()
