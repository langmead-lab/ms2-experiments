#!/usr/bin/env python3
"""Run one job of the freeze manifest: stage its inputs, run its cells, check, record, resume.

Usage:
  run_job.py JOB_ID --manifest MANIFEST --inputs-root ROOT --results-root ROOT --stage DIR
             [--platform rockfish|aws] [--no-numa] [--any-thp] [--record-traces] [--evict]

The same script runs on Rockfish and on AWS (BENCH_FREEZE_DESIGN.md D11):
- Inputs resolve by logical name under --inputs-root (a directory or an s3:// prefix) and are
  copied to --stage on local disk, each verified against its manifest sha256. Staged files are
  kept by sha256 under STAGE/inputs, so a later job on the same box reuses them; --evict frees
  files no remaining cell of this job needs when the disk would otherwise be too full.
- Cells run in an order shuffled from the manifest seed and the job id. Each cell runs its
  runs back to back (D6: three timed runs, no warmup), every command under numactl on socket 0
  (D5). A cell with a threads field (the thread-scaling matrix) runs with that many threads,
  bound with --physcpubind to its first threads physical cores, node 0's first, or with binding
  smt to both hardware threads of its first threads/2 cores, and memory on node 0. Before every
  run the page cache is dropped when running as root; otherwise each staged input of the cell is
  evicted from the cache with posix_fadvise, which needs no root but does not drop the cached
  pages other processes hold.
- With hash_output, every run's output sha256 is recorded for collect.py's cross-thread check.
- Run 1's whole output is checked against truth: ms2 MS with ms_check.cpp; ropebwt3 -l1 by the
  MS rebuilt from its SMEMs (real-hprc2 benchmark.check); ropebwt3 -lL against the SMEMs truth
  implies and an untimed --old-mem count reference (rope_validation.check_default); MONI, PHONI
  and AugPHONI with legacy_ms_validation.check_lengths; locate with offset_validation.check_smems,
  whose canonical file is kept for the cross-tool comparison of collect.py.
- After a cell's runs, its records, summary and a done marker go to the results store at once
  (D10). A rerun skips cells with a done marker, so an interrupted job loses only the cell in
  progress. Finished runs of an unfinished cell are kept on local disk and reused by a rerun on
  the same box.
- Oracle cells replay a trace given as an input. With --record-traces (Rockfish staging, D12),
  a trace still pending in the manifest is recorded by an untimed ordinary run with
  --oracle-write, checked like any run, and saved under traces/ in the results store, from where
  it is published with the indexes.

A timed-out run ends its cell: with a cell cap_s (HPRC2 phi ordinary, 3600 s) the cell is
reported as capped, otherwise as timeout. Both are final results with done markers. A failed
check or command is recorded as failed with no done marker, so a rerun tries it again.

Results: RESULTS_ROOT/MANIFEST_ID/JOB_ID/cells/CELL/{summary.json, run1.json, ..., *.stderr,
canonical.gz} and RESULTS_ROOT/MANIFEST_ID/JOB_ID/cells/CELL.done, plus job.json per attempt.
"""
import argparse
import array
import json
import os
from pathlib import Path
import platform as host_platform
import random
import shutil
import socket
import statistics
import subprocess
import sys
import time
import urllib.request

from ms2exp.lib import RESULT_FORMAT, Store, atomic_json, cell_key, file_sha256, find_job, load_manifest
from ms2exp.run import common
from ms2exp.run.timing import POLICY as TIMING_POLICY, parse as parse_phases, parse_stats
from ms2exp.validate.legacy import check_lengths
from ms2exp.validate.ms import check as check_ms_python
from ms2exp.validate.offsets import CanonicalSink, check_smems, load_truth, read_contigs, read_fasta
from ms2exp.validate.rope import check_default

QUERY_MODE = {'ms': 'min', 'psi': 'psi', 'phi': 'phi', 'dual': 'dual'}
LEGACY_OUTPUT = {'moni': 'output.lengths', 'phoni': 'reads.fa.lengths', 'augphoni': 'reads.fa.lengths.aug'}
LEGACY_PREFIX = {'moni': 'index', 'phoni': 'reference.txt', 'augphoni': 'reference.txt'}
ONE_READ_TIMEOUT = 600
# Disk kept free beyond staged inputs, for outputs, traces and references.
OUTPUT_RESERVE = 64 << 30


def log(*args):
    print(time.strftime('%H:%M:%S'), *args, flush=True)


class CellFailed(RuntimeError):
    pass


class Runner:
    def __init__(self, args):
        self.args = args
        self.manifest = load_manifest(args.manifest)
        self.job = find_job(self.manifest, args.job_id)
        self.inputs_store = Store(args.inputs_root)
        self.results = Store(args.results_root)
        self.prefix = f"{self.manifest['id']}/{self.job['id']}"
        self.stage = Path(args.stage).resolve()
        self.cache = self.stage/'inputs'
        self.work = self.stage/'work'/self.manifest['id']/self.job['id']
        self.cache.mkdir(parents=True, exist_ok=True)
        self.work.mkdir(parents=True, exist_ok=True)
        self.platform = args.platform
        self.root = os.geteuid() == 0
        self.staged = {}
        self.built = {}
        self.context = {}
        self.host = host_info()
        want = self.manifest['protocol'].get('transparent_hugepages')
        self.thp_override = False
        if want and self.job.get('pool') == 'timed' and self.host.get('thp_enabled') != want:
            if not args.any_thp:
                raise SystemExit(f"Transparent huge pages are {self.host.get('thp_enabled')!r}; the protocol "
                                 f"times every run with {want!r} (--any-thp runs anyway, as the toy may)")
            # The attempt record says the protocol's setting was overridden; summaries record the host's.
            self.thp_override = True
            log(f"transparent huge pages are {self.host.get('thp_enabled')!r}, not {want!r}; timings do not "
                'follow the protocol')

    # Selection and order.

    def cells(self):
        cells = [c for c in self.job['cells'] if self.platform in c.get('platforms', self.job['platforms'])]
        order = sorted(cells, key=lambda c: c['id'])
        random.Random(f"{self.manifest['seed']}:{self.job['id']}").shuffle(order)
        return order

    def done(self, cell):
        return self.results.exists(f"{self.prefix}/cells/{cell_key(cell['id'])}.done")

    # Staging.

    def entry(self, name):
        """An input's manifest entry; a built one takes its sha256 from its build job's outputs.json."""
        entry = self.manifest['inputs'][name]
        if entry['sha256'] != 'built':
            return entry
        if name not in self.built:
            record = f"{self.manifest['id']}/{entry['producer']}/outputs.json"
            if not self.results.exists(record):
                raise SystemExit(f"{name} is not built yet (build job {entry['producer']})")
            outputs = self.results.get_json(record)['outputs']
            self.built.update({n: dict(o, producer=entry['producer']) for n, o in outputs.items()})
        return self.built[name]

    def stage_inputs(self, cells):
        return self.stage_names({name for cell in cells for name in cell['inputs'].values()})

    def stage_names(self, names):
        needed = {name: self.entry(name) for name in names}
        pending = [n for n, e in needed.items() if e['sha256'] == 'pending']
        blocked = [n for n in pending if not (n.startswith('traces/') and self.args.record_traces)]
        if blocked:
            raise SystemExit('Inputs still pending in the manifest: ' + ', '.join(blocked))
        todo = {n: e for n, e in needed.items() if e['sha256'] != 'pending'}
        missing = sum((e.get('bytes') or 0) for e in todo.values()
                      if not (self.cache/e['sha256']/'verified').exists())
        free = shutil.disk_usage(self.cache).free
        if missing + OUTPUT_RESERVE > free:
            if self.args.evict:
                self.evict({e['sha256'] for e in todo.values()}, missing + OUTPUT_RESERVE - free)
            elif missing + OUTPUT_RESERVE > shutil.disk_usage(self.cache).free:
                raise SystemExit(f'Staging needs {missing >> 30} GiB plus reserve; {free >> 30} GiB free')
        record = {}
        for name, entry in sorted(todo.items(), key=lambda kv: kv[0]):
            started = time.time()
            path, fresh = self.stage_one(name, entry)
            self.staged[name] = path
            record[name] = dict(sha256=entry['sha256'], bytes=path.stat().st_size, fresh=fresh,
                                seconds=time.time() - started)
        # Shared libraries shipped as inputs (the pinned ropebwt3's libgomp) go on the library path.
        libs = sorted({str(self.staged[n].parent) for n in todo if n.startswith('tools/runtime/')})
        if libs:
            common.ENV['LD_LIBRARY_PATH'] = ':'.join(libs + [common.ENV.get('LD_LIBRARY_PATH', '')]).rstrip(':')
        return record

    def stage_one(self, name, entry):
        folder = self.cache/entry['sha256']
        path = folder/Path(name).name
        if (folder/'verified').exists() and path.exists():
            return path, False
        folder.mkdir(parents=True, exist_ok=True)
        tmp = folder/(path.name + '.part')
        log('stage', name)
        self.inputs_store.get(name, tmp)
        digest = file_sha256(tmp)
        if digest != entry['sha256']:
            tmp.unlink()
            raise SystemExit(f'Checksum differs for {name}: {digest} != {entry["sha256"]}')
        tmp.replace(path)
        if name.startswith('tools/'):
            path.chmod(0o755)
        (folder/'verified').write_text(json.dumps(dict(name=name, sha256=digest, at=time.time())) + '\n')
        return path, True

    def evict(self, keep, need):
        """Delete cached inputs this job does not use, least recently verified first."""
        folders = sorted((p for p in self.cache.iterdir() if p.is_dir() and p.name not in keep),
                         key=lambda p: (p/'verified').stat().st_mtime if (p/'verified').exists() else 0)
        freed = 0
        for folder in folders:
            if freed >= need:
                break
            freed += sum(f.stat().st_size for f in folder.iterdir() if f.is_file())
            shutil.rmtree(folder)
            log('evicted', folder.name)

    def file(self, cell, role):
        return self.staged[cell['inputs'][role]]

    # Shared per-read-set context.

    def workload(self, cell):
        key = cell['inputs']['offsets']
        if key not in self.context:
            bounds = array.array('Q')
            bounds.frombytes(self.file(cell, 'offsets').read_bytes())
            if sys.byteorder != 'little':
                bounds.byteswap()
            self.context[key] = dict(offsets=bounds, reads=len(bounds) - 1, bases=bounds[-1])
        return self.context[key]

    def locate_truth(self, cell):
        """Reads, MS truth and contigs for check_smems; kept for the job's other locate cells."""
        key = ('locate', cell['inputs']['reads'], cell['inputs']['lengths'])
        if key not in self.context:
            queries = read_fasta(self.file(cell, 'reads'))
            truth = load_truth(self.file(cell, 'truth'), self.file(cell, 'offsets'), queries)
            self.context[key] = dict(queries=queries, truth=truth, contigs=read_contigs(self.file(cell, 'lengths')))
        return self.context[key]

    # Commands.

    def bind(self, cell=None):
        if self.args.no_numa:
            return []
        if cell is None or 'threads' not in cell:
            return ['numactl', '--cpunodebind=0', '--membind=0']
        return ['numactl', '--physcpubind=' + ','.join(map(str, self.cpus(cell))), '--membind=0']

    def cpus(self, cell):
        """The CPUs of a threads cell: one per physical core, node 0's cores first, or with binding
        smt every hardware thread of as many cores as the thread count needs. On AWS only node 0's
        cores qualify (decision 6: socket 0 only); a Rockfish pilot node may spill onto node 1."""
        cores = [c for node, c in self.host_cores() if node == 0 or self.platform != 'aws']
        threads = cell['threads']
        if cell.get('binding') == 'smt':
            width = len(cores[0])
            if threads % width or threads // width > len(cores):
                raise CellFailed(f'{threads} threads do not fill whole cores of {width} hardware threads')
            return [c for core in cores[:threads // width] for c in core]
        if threads > len(cores):
            raise CellFailed(f'{threads} threads but {len(cores)} physical cores')
        return [core[0] for core in cores[:threads]]

    def host_cores(self):
        if 'cores' not in self.context:
            self.context['cores'] = physical_cores()
        return self.context['cores']

    def env(self, cell):
        """The command environment; a threads cell sets the OpenMP thread count to its own."""
        if 'threads' not in cell:
            return None
        return dict(common.ENV, OMP_NUM_THREADS=str(cell['threads']))

    def ms2_argv(self, cell, output, trace=None, record=None):
        argv = [self.file(cell, 'binary'), 'query', self.file(cell, 'index'), self.file(cell, 'reads'),
                '-o', output, '--threads', str(cell.get('threads', 1)), '--interleave', str(cell['interleave'])]
        if cell.get('block_bytes'):
            argv += ['--block-bytes', str(cell['block_bytes'])]
        if cell['mode'] == 'oracle':
            argv += ['--mode', 'oracle', '--oracle', trace]
        else:
            argv += ['--mode', QUERY_MODE[cell['layout']]]
            if record:
                argv += ['--oracle-write', record]
            if cell['mode'] == 'ftab':
                argv += ['--ftab', self.file(cell, 'ftab')]
        if cell['kind'] == 'locate':
            argv += ['--report', 'smem-all', '--min-smem-len', '1', '--max-smem-positions', str(cell['cap'])]
        return self.bind(cell) + [str(a) for a in argv]

    def rope_argv(self, cell, min_len=None, old_mem=False):
        options = ['-K', cell['batch']] if cell.get('batch') else []
        if old_mem:
            options += ['--old-mem']
        if cell['kind'] == 'locate':
            options += ['-l1', f"-p{cell['cap']}"]
            fmd = self.rope_dir(cell)/'index.fmd'
        else:
            length = min_len or cell['min_len']
            options += [] if (length == 19 and not old_mem) else [f'-l{length}']
            fmd = self.file(cell, 'fmd')
        return self.bind(cell) + [str(self.file(cell, 'binary')), 'mem', f"-t{cell.get('threads', 1)}", *options, str(fmd),
                              str(self.file(cell, 'reads'))]

    def rope_dir(self, cell):
        """ropebwt3 finds the SSA and lengths by the FMD's name, so locate cells get a directory of links."""
        folder = self.work/'rope'/f"s{cell['shift']}-{cell['inputs']['fmd'].replace('/', '_')}"
        if not folder.exists():
            folder.mkdir(parents=True)
            for name, role in (('index.fmd', 'fmd'), ('index.fmd.ssa', 'ssa'), ('index.fmd.len.gz', 'lengths')):
                os.symlink(self.file(cell, role), folder/name)
        return folder

    def legacy_dir(self, cell):
        folder = self.work/'legacy'/f"{cell['tool']}-{cell['reference']}"
        if not folder.exists():
            folder.mkdir(parents=True)
            for role in cell['inputs']:
                if role.startswith('bundle/'):
                    os.symlink(self.file(cell, role), folder/role.split('/', 1)[1])
        return folder

    def legacy_argv(self, cell, reads, folder):
        argv = self.bind() + [str(self.file(cell, 'binary')), str(self.legacy_dir(cell)/LEGACY_PREFIX[cell['tool']]),
                              '-p', str(reads)]
        if cell['tool'] == 'moni':
            return argv + ['-o', str(folder/'output'), '-t', '1']
        return argv + ['-g', 'naive'] + (['-t', 'compressed'] if cell['tool'] == 'augphoni' else [])

    # Cache state.

    def cold(self, cell):
        """Drop the page cache (root) or evict this cell's inputs from it (fadvise)."""
        os.sync()
        if self.root:
            Path('/proc/sys/vm/drop_caches').write_text('3\n')
            return 'drop_caches'
        if not hasattr(os, 'posix_fadvise'):
            return 'none'
        for name in set(cell['inputs'].values()):
            path = self.staged.get(name)
            if path is None:
                continue
            fd = os.open(path, os.O_RDONLY)
            try:
                os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
            finally:
                os.close(fd)
        return 'fadvise-dontneed'

    # Running.

    def execute(self, cell, argv, folder, label, timeout):
        """common.execute plus the run's cache state; returns (row, None) or (row-or-None, failure)."""
        folder.mkdir(parents=True, exist_ok=True)
        cache = self.cold(cell)
        try:
            row = common.execute(argv, folder, label, timeout=timeout, env=self.env(cell))
            row['cache'] = cache
            return row, None
        except RuntimeError as error:
            saved = folder/(label + '.command.json')
            row = json.loads(saved.read_text()) if saved.exists() else None
            if row is not None:
                row['cache'] = cache
            return row, dict(error=str(error), timed_out=bool(row and row.get('timed_out')))

    def check_ms2_output(self, cell, output):
        result = subprocess.run([str(self.file(cell, 'checker')), str(output), str(self.file(cell, 'truth')),
                                 str(self.file(cell, 'offsets'))], capture_output=True, text=True)
        lines = result.stdout.strip().splitlines()
        evidence = json.loads(lines[-1]) if lines else dict(failed=1, error=result.stderr.strip()[-2000:])
        if result.returncode or evidence.get('failed') != 0:
            raise CellFailed('ms_check: ' + json.dumps(evidence))
        evidence['validation'] = 'ms_check: every MS value equals truth'
        return evidence

    def check_output(self, cell, folder, row, output):
        """Whole-output truth check of one run; returns the evidence or raises CellFailed."""
        w = self.workload(cell)
        try:
            if cell['kind'] in ('ms', 'count') and cell['tool'].startswith('ms2'):
                return self.check_ms2_output(cell, output)
            if cell['tool'] == 'ropebwt3' and cell['kind'] == 'ms':
                if cell['min_len'] == 1:
                    evidence = check_ms_python('ropebwt3', folder, row, self.file(cell, 'truth'), w['offsets'])
                    return dict(evidence, validation='MS rebuilt from every -l1 SMEM equals truth')
                reference = self.context[('reference', cell['id'])]
                return check_default(row['stdout'], self.file(cell, 'truth'), w['offsets'],
                                     count_reference=reference, min_length=cell['min_len'])
            if cell['kind'] == 'locate':
                t = self.locate_truth(cell)
                path = self.work/'canonical'/(cell_key(cell['id']) + '.canonical.gz')
                path.parent.mkdir(parents=True, exist_ok=True)
                path.unlink(missing_ok=True)
                sink = CanonicalSink(cell['cap'], path)
                try:
                    evidence = check_smems(output, t['queries'], t['truth'], t['contigs'], cell['cap'], sink,
                                           'orbit' if cell['tool'] == 'ms2' else 'rope')
                except BaseException:
                    sink.discard()
                    raise
                evidence['canonical'] = dict(sink.close(), path=str(path))
                return evidence
        except CellFailed:
            raise
        except (ValueError, AssertionError, KeyError) as error:
            raise CellFailed(f'{type(error).__name__}: {error}')
        raise CellFailed('No check for this cell kind')

    def prepare(self, cell, cell_dir):
        """Untimed preparation: an oracle trace to record, or a ropebwt3 count reference."""
        prep = {}
        if cell['tool'].startswith('ms2') and cell['mode'] == 'oracle':
            name = cell['inputs']['trace']
            if name not in self.staged:
                prep['trace'] = self.record_trace(cell, name, cell_dir)
        if cell['tool'] == 'ropebwt3' and cell['kind'] == 'ms' and cell['min_len'] > 1:
            folder = self.work/'references'/cell_key(cell['id'])
            reference = folder/'reference.stdout'
            if not reference.exists():
                row, failure = self.execute(cell, self.rope_argv(cell, old_mem=True), folder, 'reference',
                                            cell['timeout_s'])
                if failure:
                    raise CellFailed('Count reference failed: ' + failure['error'])
                prep['reference'] = dict(argv=row['argv'], wall_seconds=row['wall_seconds'])
            self.context[('reference', cell['id'])] = reference
        return prep

    def record_trace(self, cell, name, cell_dir):
        """A pending trace: the copy another job already stored under traces/ (the timed and the
        counting cells of a read set share one), or else one recorded here and stored there."""
        trace = self.work/'traces'/name
        stored = f"{self.manifest['id']}/traces/{name}"
        if not trace.exists() and self.results.exists(stored):
            self.results.get(stored, trace)
            digest = file_sha256(trace)
            self.staged[name] = trace
            log('fetched trace', name, digest[:12])
            return dict(name=name, sha256=digest, bytes=trace.stat().st_size, stored=self.results.url(stored),
                        recorded=False)
        if not trace.exists():
            trace.parent.mkdir(parents=True, exist_ok=True)
            folder = self.work/'record'/cell_key(cell['id'])
            output = folder/'output'
            tmp = trace.with_name(trace.name + '.part')
            tmp.unlink(missing_ok=True)
            ordinary = dict(cell, mode='ordinary')
            row, failure = self.execute(cell, self.ms2_argv(ordinary, output, record=tmp), folder, 'record',
                                        cell['timeout_s'])
            if failure:
                raise CellFailed('Trace recording failed: ' + failure['error'])
            self.check_ms2_output(cell, output)
            output.unlink()
            tmp.replace(trace)
            shutil.copyfile(folder/'record.stderr', cell_dir/'record.stderr')
        digest = file_sha256(trace)
        self.results.put(trace, f"{self.manifest['id']}/traces/{name}")
        self.staged[name] = trace
        log('recorded trace', name, digest[:12])
        return dict(name=name, sha256=digest, bytes=trace.stat().st_size,
                    stored=self.results.url(f"{self.manifest['id']}/traces/{name}"))

    def run_once(self, cell, index, cell_dir):
        folder = self.work/'runs'/cell_key(cell['id'])/f'run{index}'
        if folder.exists():
            shutil.rmtree(folder)
        timeout = cell.get('cap_s') or cell['timeout_s']
        label = f'run{index}'
        record = dict(run=index, host=self.host['hostname'], started=time.time())
        check = index == 1
        if cell['tool'] in LEGACY_OUTPUT:
            return self.legacy_once(cell, index, folder, cell_dir, record)
        if cell['kind'] == 'rep':
            return self.rep_once(cell, folder, cell_dir, record)
        output = folder/'output'
        if cell['tool'].startswith('ms2'):
            trace = self.staged.get(cell['inputs'].get('trace')) if cell['mode'] == 'oracle' else None
            argv = self.ms2_argv(cell, output, trace=trace)
        else:
            argv = self.rope_argv(cell)
        row, failure = self.execute(cell, argv, folder, label, timeout)
        for suffix in ('.stderr', '.command.json'):
            if (folder/(label + suffix)).exists():
                shutil.copyfile(folder/(label + suffix), cell_dir/(label + suffix))
        if failure:
            record.update(failure, timeout_s=timeout, wall_seconds=row and row.get('wall_seconds'))
            return record
        if cell['tool'] == 'ropebwt3':
            output = Path(row['stdout'])
        stderr = Path(row['stderr']).read_text()
        record.update(argv=row['argv'], timeout_s=timeout, cache=row['cache'], wall_seconds=row['wall_seconds'],
                      user_seconds=row.get('user_seconds'), system_seconds=row.get('system_seconds'),
                      peak_rss_kib=row.get('peak_rss_kib'))
        record['phase'] = parse_phases('orbit' if cell['tool'].startswith('ms2') else 'ropebwt3', stderr)
        if cell['kind'] == 'count':
            record['ops_per_base'] = parse_stats(stderr, cell['layout'])
        if cell['kind'] == 'locate' and cell['tool'] == 'ms2':
            contigs = len(self.locate_truth(cell)['contigs'])
            if f'smem-all: positions as contig:strand:offset, {contigs} sequences' not in stderr:
                raise CellFailed('ms2 index lacks the expected sequence-ID table')
        if 'threads' in cell:
            record.update(threads=cell['threads'], cpus=self.cpus(cell) if not self.args.no_numa else None)
        if cell.get('hash_output'):
            record['output_sha256'] = file_sha256(output)
        if check:
            record['check'] = self.check_output(cell, folder, row, output)
        shutil.rmtree(folder)
        record['finished'] = time.time()
        return record

    def legacy_once(self, cell, index, folder, cell_dir, record):
        """A one-read run, then the full run, both cold; query is full wall minus one-read wall."""
        w = self.workload(cell)
        members = {}
        for member, role, timeout in (('one', 'one', ONE_READ_TIMEOUT), ('full', 'reads', cell['timeout_s'])):
            sub = folder/member
            sub.mkdir(parents=True)
            reads = sub/'reads.fa'
            os.symlink(self.file(cell, role), reads)
            row, failure = self.execute(cell, self.legacy_argv(cell, reads, sub), sub, member, timeout)
            for suffix in ('.stderr', '.command.json'):
                if (sub/(member + suffix)).exists():
                    shutil.copyfile(sub/(member + suffix), cell_dir/(f'run{index}-{member}' + suffix))
            if failure:
                record.update(failure, member=member, timeout_s=timeout, wall_seconds=row and row.get('wall_seconds'))
                return record
            out = dict(argv=row['argv'], cache=row['cache'], wall_seconds=row['wall_seconds'],
                       peak_rss_kib=row.get('peak_rss_kib'), user_seconds=row.get('user_seconds'),
                       system_seconds=row.get('system_seconds'))
            if index == 1:
                try:
                    out['check'] = check_lengths(sub/LEGACY_OUTPUT[cell['tool']], self.file(cell, 'truth'),
                                                 w['offsets'], baseline=member == 'one')
                except (ValueError, AssertionError, OSError) as error:
                    # OSError covers a tool that exits 0 without writing its output file.
                    raise CellFailed(f'{cell["tool"]} {member}: {error}')
            members[member] = out
        one, full = members['one'], members['full']
        record.update(one=one, argv=full['argv'], cache=full['cache'], timeout_s=cell['timeout_s'],
                      wall_seconds=full['wall_seconds'], user_seconds=full['user_seconds'],
                      system_seconds=full['system_seconds'],
                      peak_rss_kib=max(one.get('peak_rss_kib') or 0, full.get('peak_rss_kib') or 0) or None,
                      phase=dict(source='one-read-wall', query_definition='full run wall minus one-read run wall',
                                 load_seconds=one['wall_seconds'],
                                 query_seconds=full['wall_seconds'] - one['wall_seconds']))
        if index == 1:
            record['check'] = dict(full['check'], one_read=one['check'])
        shutil.rmtree(folder)
        record['finished'] = time.time()
        return record

    def rep_once(self, cell, folder, cell_dir, record):
        """A1: fig_rep.py run (psi, phi and dual with TMS_REP_LOG), which fails unless the modes agree."""
        dest = folder/'rep'
        argv = [*self.bind(), sys.executable, '-m', 'ms2exp.figures.fig_rep', 'run', '--binary', str(self.file(cell, 'binary')),
                '--reads', str(self.file(cell, 'reads')), '--destination', str(dest),
                '--interleave', str(cell.get('interleave', 32))]
        for mode in ('psi', 'phi', 'dual'):
            argv += [f'--index-{mode}', str(self.file(cell, 'index_' + mode))]
        row, failure = self.execute(cell, argv, folder, 'rep', cell['timeout_s'])
        for suffix in ('.stderr', '.command.json'):
            if (folder/('rep' + suffix)).exists():
                shutil.copyfile(folder/('rep' + suffix), cell_dir/('rep' + suffix))
        if failure:
            record.update(failure, timeout_s=cell['timeout_s'])
            return record
        rep = json.loads((dest/'summary.json').read_text())
        times = rep['times']
        # The replogs and run.json are what fig_rep_violin.py and fig_rep_tikz.py read.
        shutil.copytree(dest, cell_dir/'rep', dirs_exist_ok=True)
        record.update(argv=row['argv'], cache=row['cache'], wall_seconds=row['wall_seconds'],
                      peak_rss_kib=row.get('peak_rss_kib'), timeout_s=cell['timeout_s'],
                      check=dict(outputs_identical=rep['outputs_identical'], aligned=rep['aligned'],
                                 log_sums_match_stats=True, bases=rep['bases'], repositions=rep['repositions']),
                      phase=dict(source='fig_rep.py', query_definition='sum of the three modes\' query_s',
                                 query_seconds=sum(t['query_s'] for t in times.values()),
                                 load_seconds=sum(t['index_load_s'] for t in times.values())),
                      rep=dict(repositions_per_base=rep['repositions_per_base'], times=times))
        shutil.rmtree(folder)
        record['finished'] = time.time()
        return record

    def run_cell(self, cell):
        key = cell_key(cell['id'])
        cell_dir = self.work/'results'/key
        cell_dir.mkdir(parents=True, exist_ok=True)
        runs = cell.get('runs', self.manifest['protocol']['runs'])
        summary = self.base_summary(cell)
        try:
            summary['preparation'] = self.prepare(cell, cell_dir)
            records = []
            for index in range(1, runs + 1):
                saved = cell_dir/f'run{index}.json'
                if saved.exists():
                    records.append(json.loads(saved.read_text()))
                    log(cell['id'], f'run{index} reused')
                    continue
                log(cell['id'], f'run{index}')
                record = self.run_once(cell, index, cell_dir)
                atomic_json(saved, record)
                records.append(record)
                if 'error' in record:
                    break
            summary.update(self.summarize(cell, records))
        except CellFailed as error:
            summary.update(status='failed', error=str(error))
        summary['finished'] = time.time()
        canonical = self.work/'canonical'/(key + '.canonical.gz')
        if canonical.exists():
            shutil.copyfile(canonical, cell_dir/'canonical.gz')
        atomic_json(cell_dir/'summary.json', summary)
        self.results.put_tree(cell_dir, f'{self.prefix}/cells/{key}')
        if summary['status'] in ('checked', 'capped', 'timeout'):
            self.results.put_json(dict(cell=cell['id'], status=summary['status'], finished=summary['finished'],
                                       summary_sha256=file_sha256(cell_dir/'summary.json')),
                                  f'{self.prefix}/cells/{key}.done')
        log(cell['id'], summary['status'], summary.get('median_query_seconds'), summary.get('error', ''))
        return summary

    def base_summary(self, cell):
        # A1 read sets have no offsets file; fig_rep.py counts their bases.
        w = self.workload(cell) if 'offsets' in cell['inputs'] else dict(reads=None, bases=None)
        fig_id = cell['id'].split('/', 1)[1]
        out = dict(format=RESULT_FORMAT, id=fig_id, cell=cell['id'], group=cell['group'],
                   figures=cell['figures'], tool=cell['tool'], kind=cell['kind'],
                   manifest_id=self.manifest['id'], job=self.job['id'], platform=self.platform,
                   host=self.host, timing_policy=TIMING_POLICY, protocol=self.manifest['protocol'],
                   ms2=self.manifest['ms2'], ropebwt3=self.manifest['ropebwt3'],
                   reads=w['reads'], bases=w['bases'], reference=cell['reference'], read_set=cell['read_set'],
                   inputs={role: dict(name=name, sha256=self.entry(name)['sha256'])
                           for role, name in cell['inputs'].items()},
                   numa=not self.args.no_numa, root=self.root, started=time.time())
        for field in ('layout', 'mode', 'interleave', 'block_bytes', 'cap', 'min_len', 'shift', 'batch', 'cap_s',
                      'threads', 'binding'):
            if field in cell:
                out[field] = cell[field]
        return out

    def summarize(self, cell, records):
        failed = [r for r in records if 'error' in r]
        if failed:
            f = failed[0]
            if f.get('timed_out'):
                status = 'capped' if cell.get('cap_s') else 'timeout'
                return dict(status=status, runs=records, timeout_s=f.get('timeout_s'),
                            note=f"run {f['run']} reached its {f.get('timeout_s')} s limit")
            raise CellFailed(f"run {f['run']}: {f['error']}")
        query = [r['phase']['query_seconds'] for r in records]
        load = [r['phase'].get('load_seconds') for r in records]
        med = statistics.median(query)
        if cell['kind'] == 'rep':
            return dict(status='checked', runs=records, check=records[0]['check'], rep=records[0]['rep'],
                        bases=records[0]['check']['bases'], query_seconds=query, median_query_seconds=med,
                        load_seconds=load, median_load_seconds=statistics.median(load),
                        peak_rss_kib=records[0].get('peak_rss_kib'), rep_dir='rep')
        bases = self.workload(cell)['bases']
        out = dict(status='checked', runs=records, check=records[0].get('check'),
                   query_seconds=query, median_query_seconds=med,
                   load_seconds=load, median_load_seconds=statistics.median(load) if None not in load else None,
                   wall_seconds=[r['wall_seconds'] for r in records],
                   peak_rss_kib=max((r.get('peak_rss_kib') or 0) for r in records) or None,
                   query_spread_fraction=(max(query) - min(query)) / med if med > 0 else None,
                   cache=sorted({r.get('cache') for r in records}))
        if med > 0:
            out.update(mbases_per_second=bases / med / 1e6, ns_per_base=med / bases * 1e9)
        if cell['kind'] == 'locate':
            # As explore_offsets.cell_summary: locations are the reported positions, checked in bounds.
            locations = out['check'].get('reported_locations_checked')
            out.update(locations=locations, smems=out['check'].get('intervals'))
            if locations is not None and med > 0:
                out['locations_per_second'] = locations / med
        if cell['kind'] == 'count':
            out['ops_per_base'] = records[0].get('ops_per_base')
        if cell.get('hash_output'):
            out['output_sha256'] = sorted({r['output_sha256'] for r in records})
        return out

    def run(self):
        attempt = dict(job=self.job['id'], manifest_id=self.manifest['id'], host=self.host, platform=self.platform,
                       started=time.time(), argv=sys.argv, root=self.root)
        if self.thp_override:
            attempt['thp_override'] = True
        cells = self.cells()
        todo = [c for c in cells if not self.done(c)]
        attempt.update(cells=len(cells), already_done=len(cells) - len(todo),
                       order=[c['id'] for c in cells])
        log(f"{self.job['id']}: {len(todo)} of {len(cells)} cells to run")
        if todo:
            attempt['staging'] = self.stage_inputs(todo)
        results = {}
        for cell in todo:
            results[cell['id']] = self.run_cell(cell)['status']
        attempt.update(finished=time.time(), results=results)
        stamp = time.strftime('%Y%m%dT%H%M%S')
        self.results.put_json(attempt, f"{self.prefix}/attempts/{stamp}-{self.host['hostname']}.json")
        remaining = [c['id'] for c in cells if not self.done(c)]
        if not remaining:
            self.results.put_json(dict(job=self.job['id'], finished=time.time(), cells=len(cells)),
                                  f'{self.prefix}/job.done')
        if not self.args.keep_work and not remaining:
            shutil.rmtree(self.work, ignore_errors=True)
        log(f"{self.job['id']}: {len(remaining)} cells without a done marker")
        return 0 if not remaining else 1


def host_info():
    info = dict(hostname=socket.gethostname(), machine=host_platform.machine(), kernel=host_platform.release())
    try:
        cpu = next(l for l in Path('/proc/cpuinfo').read_text().splitlines() if l.startswith('model name'))
        info['cpu'] = cpu.split(':', 1)[1].strip()
    except (OSError, StopIteration):
        pass
    for key, path in (('instance_type', 'instance-type'), ('instance_id', 'instance-id'),
                      ('availability_zone', 'placement/availability-zone'), ('lifecycle', 'instance-life-cycle')):
        value = imds(path)
        if value:
            info[key] = value
    if os.environ.get('SLURM_JOB_ID'):
        info['slurm_job_id'] = os.environ['SLURM_JOB_ID']
    # The bracketed choice in each sysfs file, for example 'always' from '[always] madvise never'.
    for key in ('enabled', 'defrag'):
        try:
            text = Path(f'/sys/kernel/mm/transparent_hugepage/{key}').read_text()
            info[f'thp_{key}'] = text.split('[', 1)[1].split(']', 1)[0]
        except (OSError, IndexError):
            pass
    return info


def cpu_list(text):
    """CPU numbers of a sysfs list such as '0-3,8,10-11'."""
    out = []
    for part in text.strip().split(','):
        if part:
            a, _, b = part.partition('-')
            out += range(int(a), int(b or a) + 1)
    return out


def physical_cores(sys_root='/sys/devices/system'):
    """(node, hardware threads) of each physical core, node 0's cores first, then node 1's and so
    on, each node's ordered by its first hardware thread."""
    root = Path(sys_root)
    nodes = sorted((int(p.name[4:]), p) for p in (root/'node').glob('node[0-9]*'))
    out = []
    for number, node in nodes:
        cores = {}
        for cpu in cpu_list((node/'cpulist').read_text()):
            siblings = cpu_list((root/'cpu'/f'cpu{cpu}'/'topology'/'thread_siblings_list').read_text())
            cores.setdefault(tuple(sorted(siblings)), sorted(siblings))
        out += [(number, c) for c in sorted(cores.values())]
    return out


def imds(path):
    """One EC2 instance-metadata value (IMDSv2), or None off EC2."""
    try:
        token = urllib.request.urlopen(urllib.request.Request(
            'http://169.254.169.254/latest/api/token', method='PUT',
            headers={'X-aws-ec2-metadata-token-ttl-seconds': '60'}), timeout=0.5).read().decode()
        return urllib.request.urlopen(urllib.request.Request(
            'http://169.254.169.254/latest/meta-data/' + path,
            headers={'X-aws-ec2-metadata-token': token}), timeout=0.5).read().decode()
    except OSError:
        return None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('job_id')
    p.add_argument('--manifest', required=True)
    p.add_argument('--inputs-root', required=True)
    p.add_argument('--results-root', required=True)
    p.add_argument('--stage', required=True, help='local disk for staged inputs and work files')
    p.add_argument('--platform', choices=('rockfish', 'aws'), required=True)
    p.add_argument('--no-numa', action='store_true', help='run without numactl (tests, or hosts without it)')
    p.add_argument('--any-thp', action='store_true',
                   help='run timed jobs under any transparent huge page setting, recording the override')
    p.add_argument('--record-traces', action='store_true', help='record oracle traces still pending (D12)')
    p.add_argument('--evict', action='store_true', help='free cached inputs other jobs staged when short of disk')
    p.add_argument('--keep-work', action='store_true')
    sys.exit(Runner(p.parse_args()).run())


if __name__ == '__main__':
    main()
