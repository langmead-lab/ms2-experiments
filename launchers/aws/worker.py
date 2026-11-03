#!/usr/bin/env python3
"""AWS box loop of the freeze benchmark: claim a job from the S3 queue, run it, repeat, then leave.

Usage (userdata.sh runs it as root on each box):
  worker.py --manifest M --root s3://bucket/prefix --work /work [--terminate-when-empty]

The queue is the manifest itself: every job whose platforms include aws and that has no
job.done marker under ROOT/results/MANIFEST_ID/JOB/. A box claims a job by creating
JOB/claim.json with an S3 conditional write (If-None-Match), so two boxes never run one job.
While a job runs, the box rewrites JOB/heartbeat.json every few minutes; a claim whose heartbeat
is older than STALE_S is taken over (its box was reclaimed without notice). On a spot
interruption notice the box releases its claim and stops its job at once, so the box the Auto
Scaling group launches in its place resumes the job from its per-cell done markers (D10).

All-AWS manifests (ms2exp.manifest --all-aws) add build jobs, which run ms2exp.build.job and publish
indexes and binaries that other jobs use (inputs with sha256 "built"). A job is offered only when
every build it depends on has its job.done, so the queue releases work as soon as each index
exists. Build jobs are offered first. After each job, and while waiting, a box raises the group's
desired capacity to the number of jobs running plus the jobs ready and unclaimed (never above the
group's maximum), so new boxes start as soon as work exists for them.

When nothing is ready, a box waits while some unfinished job may still become ready, up to
--max-idle-minutes, then uploads its log and, with --terminate-when-empty, removes itself from
the Auto Scaling group with a decrement of the desired capacity, so the group shrinks instead of
replacing it. Results and the staged-input cache live on /work.
"""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time

from ms2exp.lib import Store, load_manifest
from ms2exp.run.runner import imds

HEARTBEAT_S = 180
STALE_S = 1200
POLL_S = 60


def log(*args):
    print(time.strftime('%Y-%m-%dT%H:%M:%S'), *args, flush=True)


class Worker:
    def __init__(self, args):
        self.args = args
        self.manifest = load_manifest(args.manifest)
        self.root = args.root.rstrip('/')
        self.results = Store(self.root + '/results')
        self.instance = imds('instance-id') or os.uname().nodename
        self.instance_type = imds('instance-type')
        self.current = None
        self.process = None
        self.interrupted = threading.Event()
        self.stop = threading.Event()

    def key(self, job, name):
        return f"{self.manifest['id']}/{job['id']}/{name}"

    def jobs(self):
        """AWS jobs, build jobs first."""
        jobs = [j for j in self.manifest['jobs'] if 'aws' in j['platforms']]
        return sorted(jobs, key=lambda j: j.get('kind') != 'build')

    def ready(self, job, done):
        """Whether every build this job's inputs come from has finished (done: set of finished job ids)."""
        for name in job['inputs']:
            entry = self.manifest['inputs'][name]
            if entry['sha256'] == 'pending' and not name.startswith('traces/'):
                return False
            if entry['sha256'] == 'built' and entry['producer'] not in done:
                return False
        return True

    def snapshot(self):
        """Finished, claimed and ready-unclaimed job ids, from the results store."""
        done, claimed, ready = set(), set(), []
        present = set(self.results.names(self.manifest['id']))
        for job in self.jobs():
            if self.key(job, 'job.done') in present:
                done.add(job['id'])
            elif self.key(job, 'claim.json') in present and not self.stale(job):
                claimed.add(job['id'])
        for job in self.jobs():
            if job['id'] not in done and job['id'] not in claimed and self.ready(job, done):
                ready.append(job['id'])
        return done, claimed, ready

    def scale_up(self, claimed, ready):
        """Raise the group's desired capacity to running plus ready jobs (instance-role credentials)."""
        group = self.group()
        if group is None:
            return
        want = min(group['MaxSize'], len(claimed) + len(ready))
        if want > group['DesiredCapacity']:
            log('desired capacity', group['DesiredCapacity'], '->', want)
            autoscaling('set-desired-capacity', '--auto-scaling-group-name', group['AutoScalingGroupName'],
                        '--desired-capacity', str(want), '--honor-cooldown')

    def group(self):
        if not hasattr(self, '_group_name'):
            found = autoscaling('describe-auto-scaling-instances', '--instance-ids', self.instance)
            rows = (found or {}).get('AutoScalingInstances', [])
            self._group_name = rows[0]['AutoScalingGroupName'] if rows else None
        if self._group_name is None:
            return None
        found = autoscaling('describe-auto-scaling-groups', '--auto-scaling-group-names', self._group_name)
        groups = (found or {}).get('AutoScalingGroups', [])
        return groups[0] if groups else None

    def stale(self, job):
        """Whether job's claim belongs to another box that has written neither the claim nor a
        heartbeat for STALE_S (it was reclaimed without notice). A heartbeat left by an earlier
        holder does not count for the current one."""
        try:
            current = self.results.get_json(self.key(job, 'claim.json'))
        except RuntimeError:
            return False
        holder, last = current.get('instance'), current.get('claimed', 0)
        try:
            beat = self.results.get_json(self.key(job, 'heartbeat.json'))
            if beat.get('instance') == holder:
                last = max(last, beat.get('time', 0))
        except RuntimeError:
            pass
        return holder != self.instance and time.time() - last >= STALE_S

    def claim(self, job):
        if self.results.exists(self.key(job, 'job.done')):
            return False
        claim = dict(instance=self.instance, instance_type=self.instance_type, claimed=time.time())
        if self.results.put_if_absent(claim, self.key(job, 'claim.json')):
            return True
        if not self.stale(job):
            return False
        log('taking over stale claim of', job['id'])
        self.results.delete(self.key(job, 'claim.json'))
        return self.results.put_if_absent(claim, self.key(job, 'claim.json'))

    def release(self, job):
        self.results.delete(self.key(job, 'claim.json'))

    def heartbeat(self):
        while not self.stop.wait(HEARTBEAT_S):
            job = self.current
            if job is not None:
                try:
                    self.results.put_json(dict(instance=self.instance, time=time.time(), job=job['id']),
                                          self.key(job, 'heartbeat.json'))
                except RuntimeError as error:
                    log('heartbeat failed:', error)

    def watch_spot(self):
        """Poll for a spot interruption notice (two minutes of warning)."""
        while not self.stop.wait(5):
            if imds('spot/instance-action'):
                log('spot interruption notice')
                self.interrupted.set()
                job = self.current
                if job is not None:
                    self.release(job)
                if self.process is not None:
                    os.killpg(self.process.pid, signal.SIGTERM)
                return

    def run_job(self, job):
        module = 'ms2exp.build.job' if job.get('kind') == 'build' else 'ms2exp.run.runner'
        argv = [sys.executable, '-m', module, job['id'], '--manifest', self.args.manifest,
                '--inputs-root', self.args.inputs_root or self.root + '/inputs', '--results-root', self.root + '/results',
                '--stage', str(Path(self.args.work)/'stage'), '--platform', 'aws', '--evict']
        if module == 'ms2exp.run.runner':
            # Without a Rockfish dry run (all-AWS), oracle cells record their own traces untimed first.
            argv.append('--record-traces')
        log('run', job['id'])
        self.process = subprocess.Popen(argv, start_new_session=True)
        status = self.process.wait()
        self.process = None
        log('job', job['id'], 'exit', status)
        return status

    def main(self):
        threading.Thread(target=self.heartbeat, daemon=True).start()
        threading.Thread(target=self.watch_spot, daemon=True).start()
        failures = {}
        idle_since = None
        while not self.interrupted.is_set():
            done, claimed, ready = self.snapshot()
            self.scale_up(claimed, ready)
            offered = [j for j in self.jobs() if j['id'] in ready and failures.get(j['id'], 0) < 2]
            job = next((j for j in offered if self.claim(j)), None)
            if job is None:
                waiting = [j for j in self.jobs() if j['id'] not in done and failures.get(j['id'], 0) < 2]
                idle_since = idle_since or time.time()
                if not waiting or time.time() - idle_since > self.args.max_idle_minutes * 60:
                    break
                log('nothing ready;', len(waiting), 'jobs unfinished; waiting')
                if self.interrupted.wait(POLL_S):
                    break
                continue
            idle_since = None
            self.current = job
            status = self.run_job(job)
            self.current = None
            if self.interrupted.is_set():
                break
            if status != 0:
                # Cells without a done marker failed; one more attempt here, then leave it to another box.
                failures[job['id']] = failures.get(job['id'], 0) + 1
            self.release(job)
        self.stop.set()
        log_path = Path(self.args.work)/'worker.out'
        if log_path.exists():
            self.results.put(log_path, f"{self.manifest['id']}/boxes/{self.instance}.log")
        if self.interrupted.is_set():
            return
        log('queue empty')
        if self.args.terminate_when_empty:
            terminate(self.instance)


def autoscaling(*args):
    """One autoscaling CLI call with the instance role (JSON result, or None on failure)."""
    env = {k: v for k, v in os.environ.items() if k != 'AWS_PROFILE'}
    result = subprocess.run(['aws', 'autoscaling', *args, '--output', 'json'], env=env, capture_output=True,
                            text=True)
    if result.returncode:
        log('autoscaling', args[0], 'failed:', result.stderr.strip())
        return None
    return json.loads(result.stdout) if result.stdout.strip() else {}


def terminate(instance):
    """Leave the Auto Scaling group with a decrement (instance-role credentials); else power off."""
    env = {k: v for k, v in os.environ.items() if k != 'AWS_PROFILE'}
    result = subprocess.run(['aws', 'autoscaling', 'terminate-instance-in-auto-scaling-group', '--instance-id',
                             instance, '--should-decrement-desired-capacity'], env=env, capture_output=True, text=True)
    log('terminate:', result.returncode, result.stdout.strip(), result.stderr.strip())
    if result.returncode:
        subprocess.run(['shutdown', '-h', 'now'])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest', required=True)
    p.add_argument('--root', required=True, help='bench S3 root: inputs/, results/ and code/ below it')
    p.add_argument('--inputs-root', help='S3 prefix of the inputs (default ROOT/inputs)')
    p.add_argument('--work', required=True)
    p.add_argument('--terminate-when-empty', action='store_true')
    p.add_argument('--max-idle-minutes', type=float, default=20,
                   help='how long to wait for a job to become ready before leaving')
    Worker(p.parse_args()).main()


if __name__ == '__main__':
    main()
