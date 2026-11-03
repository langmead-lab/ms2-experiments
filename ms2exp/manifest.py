#!/usr/bin/env python3
"""Generate the freeze benchmark job manifest (BENCH_FREEZE_DESIGN.md D8, D9, D11).

The manifest is the committed, declarative list of jobs. Each job is one index group and
one read set, and lists its cells. Every input is a logical name with its sha256 and size,
taken from an inputs catalog (ms2exp.catalog hash). The ms2 source (repo URL plus commit), the
ropebwt3 pin and the D6 protocol are recorded once at the top and in every result.

Matrices:
  freeze         the paper's cells: Fig 3 timing and counting, Fig 4/A2, Fig 5/A5, Fig 6/A3/A4/A6
  toy            mtb152-sized cells that exercise every code path of the runner (toy_prep.sh)
  threads        thread scaling (THREAD_SCALING_PLAN.md): HPRC2 short and ONT 1 Gbase sets, ms2
                 min, psi, phi and dual and ropebwt3 -l1, at 1 to 32 threads on physical cores and
                 64 on both hardware threads of 32 cores; AWS only
  threads-pilot  the threads cells at 1 and 32 threads on the Fig 5 60 Mb sets, which are prefixes
                 of the 1 Gbase sets, and the 64-thread SMT cell (decision 10): ms2 on one AWS box,
                 ropebwt3 on Rockfish

Usage:
  python -m ms2exp.manifest --matrix freeze --inputs inputs.json --ms2-repo URL --ms2-commit SHA \
      [--allow-pending] [--seed N] -o manifest.json

An input missing from the catalog is an error unless --allow-pending, which records it with
sha256 "pending" (useful before the final ms2 indexes exist); the runner refuses to stage a
pending input, except oracle traces, which a Rockfish run records itself (D12).

Cell ids are GROUP/FIG_ID. GROUP names a summary that ms2exp.figures.paper_figures reads (ms2exp.collect
writes GROUP/summary.json) and FIG_ID is the id paper_figures looks up in it, as in the
exploratory runs: {layout}-{mode} for ms2 MS, rope-l1, rope-default (length 19) and rope-lL,
moni, phoni, augphoni, and cap{N}-{layout}-{ordinary|ftab} and cap{N}-rope-s{shift} for
locate. Layout ms is ms2's min layout (ms2 query --mode min).
"""
import argparse
import json
from pathlib import Path

from ms2exp.lib import MANIFEST_FORMAT, json_sha256, now

LAYOUTS = ('ms', 'psi', 'phi', 'dual')
QUERY_MODE = {'ms': 'min', 'psi': 'psi', 'phi': 'phi', 'dual': 'dual'}
ROPE_LENGTHS = (1, 19, 25, 30, 40)
CAPS = (1, 10, 100)
HPRC2_MS_SETS = ('human_short', 'human_pacbio', 'human_ont_modern', 'guccione_skin')
# Locate read sets (D9, D13.5): whole-read prefixes of the 60 Mb sets, short 6 Mb, PacBio 20 Mb,
# ONT 10 Mb, skin 0.5 Mb (d13_prep.py cap).
HPRC2_LOCATE_SETS = {'human_short': '6mb', 'human_pacbio': '20mb', 'human_ont_modern': '10mb',
                     'guccione_skin': '500kb'}
CHR22_SETS = ('simshort', 'simpacbio', 'simont')
PROTOCOL = dict(
    runs=3, warmup_runs=0, report='median',
    drop_page_cache='before every run, when root (echo 3 > /proc/sys/vm/drop_caches)',
    numa=dict(cpunodebind=0, membind=0), threads=1,
    # ms2 does not request huge pages itself, so its large indexes sit on 4 KB pages unless THP is
    # always on. A timed job refuses to start under any other setting.
    transparent_hugepages='always',
    truth_check='run 1 of every cell, whole output',
    order='cells in an order shuffled from the manifest seed and the job id; each cell runs its three runs '
          'back to back, so a done cell is final (D10)',
    time_basis=dict(ms2='phase timer: query = total_s (read, parse, query, format, write), load = '
                        'index_load_s + ftab_s',
                    ropebwt3='log stamps: load = last loaded stamp, query = last processed stamp minus load',
                    legacy='whole-process wall of the full run minus the one-read run, both cold'),
    design='BENCH_FREEZE_DESIGN.md D5, D6, D10')
# Thread scaling (THREAD_SCALING_PLAN.md decisions 5 to 7): physical cores of socket 0, then both
# hardware threads of all 32 of them; ms2 splits the reads into 2 MiB blocks at every thread
# count, and ropebwt3 keeps its default batch.
THREADS = (1, 2, 4, 8, 16, 32, 64)
PILOT_THREADS = (1, 32)
THREAD_SETS = ('human_short', 'human_ont_modern')
THREAD_BLOCK = 2 << 20
THREADS_PROTOCOL = dict(
    PROTOCOL, threads='per cell',
    numa=dict(physcpubind='the first T physical cores, node 0 first; T = 64 on both hardware threads of the '
                          'first 32 cores (binding smt)', membind=0),
    cross_thread_check='every run of every cell records its output sha256, which must be the same at every '
                       'thread count of a tool, mode and read set (ms2exp.collect)',
    design='BENCH_FREEZE_DESIGN.md D5, D6, D10; THREAD_SCALING_PLAN.md')
# Per-run limits in seconds. HPRC2 phi ordinary keeps the 3600 s cap and is reported as capped.
TIMEOUT = dict(ms2=7200, ropebwt3=7200, legacy=7200, count=14400, locate=14400)
PHI_CAP_S = 3600


class Builder:
    def __init__(self, catalog, allow_pending, ms2_commit, rope_tag, all_aws=False):
        self.catalog, self.allow_pending = catalog, allow_pending
        self.ms2_commit, self.rope_tag = ms2_commit, rope_tag
        self.inputs, self.jobs, self.missing = {}, [], []
        # Logical name -> id of the build job that produces it (all-AWS mode only).
        self.produced = {}
        self.all_aws = all_aws
        # Counting, A1 and D9 cells run on Rockfish in the staged plan (D7, D9); in all-AWS mode
        # they run on the AWS boxes too.
        self.untimed_platforms = ['rockfish', 'aws'] if all_aws else ['rockfish']
        self.d9_platforms = self.untimed_platforms

    def input(self, name, pending_ok=False):
        if name not in self.inputs and name in self.produced:
            # Its checksum exists only once the build ran: the runner reads it from the build job's
            # outputs.json, and worker.py offers the job only after that build is done.
            self.inputs[name] = dict(sha256='built', bytes=None, producer=self.produced[name])
        if name not in self.inputs:
            entry = self.catalog.get(name)
            if entry is None:
                if not (self.allow_pending or pending_ok):
                    self.missing.append(name)
                entry = dict(sha256='pending', bytes=None)
            self.inputs[name] = dict(sha256=entry['sha256'], bytes=entry.get('bytes'))
        return name

    def job(self, job_id, cells, platforms=('rockfish', 'aws'), note=None, pool='timed'):
        ids = [c['id'] for c in cells]
        if len(set(ids)) != len(ids):
            raise ValueError('Duplicate cell id in ' + job_id)
        inputs = sorted({v for c in cells for v in c['inputs'].values()})
        job = dict(id=job_id, kind='cells', pool=pool, platforms=list(platforms), inputs=inputs, cells=cells)
        if note:
            job['note'] = note
        self.jobs.append(job)

    # Tools.
    def ms2(self, stats=False, replog=False):
        name = 'ms2-replog' if replog else ('ms2-stats' if stats else 'ms2-bench')
        return self.input(f'tools/ms2-{self.ms2_commit[:7]}/{name}')

    def rope(self):
        return self.input(f'tools/ropebwt3-{self.rope_tag}/ropebwt3')

    def runtime(self, soname='libgomp.so.1'):
        """A shared library shipped with the pinned binaries (ms2exp.build.tools), put on LD_LIBRARY_PATH:
        by default the OpenMP runtime the pinned ropebwt3 links."""
        return self.input(f'tools/runtime/{soname}')

    def checker(self):
        return self.input('tools/ms_check')

    def workload(self, ref, rset, one=False):
        files = dict(reads=self.input(f'reads/{ref}/{rset}/reads.fa'),
                     truth=self.input(f'reads/{ref}/{rset}/truth.u32'),
                     offsets=self.input(f'reads/{ref}/{rset}/offsets.u64'))
        if one:
            files['one'] = self.input(f'reads/{ref}/{rset}/one.fa')
        return files

    # Cells.
    def ms2_ms(self, group, ref, rset, layout, mode, index_name=None, interleave=32, block_bytes=None,
               figures=(), cap_s=None, fig_id=None, platforms=None, timeout=None):
        """ms2 matching statistics; mode ordinary, oracle or ftab."""
        index = self.input(f'ms2/{ref}/{index_name or default_index(layout, ref)}')
        inputs = dict(binary=self.ms2(), checker=self.checker(), index=index, **self.workload(ref, rset))
        if mode == 'ftab':
            inputs['ftab'] = self.input(index + '.ftab')
        if mode == 'oracle':
            inputs['trace'] = self.input(f'traces/{ref}/{rset}/{layout}.oracle', pending_ok=True)
        cell = dict(id=f'{group}/{fig_id or layout + "-" + mode}', group=group, figures=list(figures),
                    tool='ms2', kind='ms', layout=layout, mode=mode, interleave=interleave,
                    reference=ref, read_set=rset, inputs=inputs, timeout_s=timeout or TIMEOUT['ms2'])
        if block_bytes:
            cell['block_bytes'] = block_bytes
        if cap_s:
            cell['cap_s'] = cap_s
        if platforms:
            cell['platforms'] = list(platforms)
        return cell

    def ms2_count(self, group, ref, rset, layout, mode, interleave=32, block_bytes=None, figures=()):
        """Counting build (MS_STATS, TMS_STATS): one untimed run whose stats line gives operations per base."""
        cell = self.ms2_ms(group, ref, rset, layout, mode, interleave=interleave, block_bytes=block_bytes,
                           figures=figures, timeout=TIMEOUT['count'])
        cell['inputs']['binary'] = self.ms2(stats=True)
        cell.update(tool='ms2-stats', kind='count', runs=1, platforms=self.untimed_platforms)
        return cell

    def ms2_rep(self, group, ref, rset, fig_id, reads=None, figures=('A1',)):
        """A1: fig_rep.py run with the TMS_REP_LOG build, psi, phi and dual on one read set (one untimed run).

        fig_rep.py checks its own output: the three modes' MS must be identical and their logs must agree
        record for record and with the stats line. No MS truth exists for these read sets."""
        index = {m: self.input(f'ms2/{ref}/{default_index(m, ref)}') for m in ('psi', 'phi', 'dual')}
        inputs = dict(binary=self.ms2(replog=True), reads=self.input(reads or f'reads/{ref}/{rset}/reads.fa'),
                      index_psi=index['psi'], index_phi=index['phi'], index_dual=index['dual'])
        return dict(id=f'{group}/{fig_id}', group=group, figures=list(figures), tool='ms2-replog', kind='rep',
                    runs=1, interleave=32, reference=ref, read_set=rset, inputs=inputs,
                    timeout_s=TIMEOUT['count'], platforms=self.untimed_platforms)

    # Build jobs (all-AWS mode).
    def build(self, job_id, steps, outputs, note=None):
        """A job that runs build commands and publishes its outputs as inputs (ms2exp.build.job).

        Step argv tokens: {in:NAME} a staged input, {out:NAME} the output being written (also usable as
        the input of a later step), {outroot} the directory outputs are written under by logical name,
        {threads} the box's CPU count, {code} the freeze code directory, {work} a scratch directory."""
        names = {t[4:-1] for step in steps for t in step['argv'] if t.startswith('{in:') and t.endswith('}')}
        for name in outputs:
            if name in self.produced:
                raise ValueError(f'{name} is produced twice')
            self.produced[name] = job_id
        inputs = sorted(self.input(n) for n in names)
        job = dict(id=job_id, kind='build', pool='build', platforms=['aws'], inputs=inputs,
                   outputs=sorted(outputs), steps=steps, cells=[])
        if note:
            job['note'] = note
        self.jobs.append(job)

    def rope_ms(self, group, ref, rset, length, batch=None, figures=()):
        fig_id = 'rope-l1' if length == 1 else ('rope-default' if length == 19 else f'rope-l{length}')
        inputs = dict(binary=self.rope(), runtime=self.runtime(), fmd=self.input(f'ropebwt3/{ref}/index.fmd'),
                      **self.workload(ref, rset))
        cell = dict(id=f'{group}/{fig_id}', group=group, figures=list(figures), tool='ropebwt3', kind='ms',
                    min_len=length, reference=ref, read_set=rset, inputs=inputs, timeout_s=TIMEOUT['ropebwt3'])
        if batch:
            cell['batch'] = batch
        return cell

    def legacy(self, group, ref, rset, tool, figures=()):
        files = LEGACY_FILES[tool]
        inputs = dict(binary=self.input(f'tools/{tool}/{LEGACY_BINARY[tool]}'), **self.workload(ref, rset, one=True))
        for name in files:
            inputs['bundle/' + name] = self.input(f'{tool}/{ref}/{name}')
        for soname in LEGACY_RUNTIME:
            inputs['runtime/' + soname] = self.runtime(soname)
        return dict(id=f'{group}/{tool}', group=group, figures=list(figures), tool=tool, kind='ms',
                    reference=ref, read_set=rset, inputs=inputs, timeout_s=TIMEOUT['legacy'])

    def ms2_locate(self, group, ref, rset, layout, cap, mode='ordinary', index_name=None, figures=(),
                   interleave=32):
        index = self.input(f'ms2/{ref}/{index_name or locate_index(layout)}')
        inputs = dict(binary=self.ms2(), index=index, lengths=self.input(f'ropebwt3/{ref}/index.fmd.len.gz'),
                      **self.workload(ref, rset))
        if mode == 'ftab':
            inputs['ftab'] = self.input(index + '.ftab')
        return dict(id=f'{group}/cap{cap}-{layout}-{mode}', group=group, figures=list(figures), tool='ms2',
                    kind='locate', layout=layout, mode=mode, cap=cap, interleave=interleave, reference=ref,
                    read_set=rset, inputs=inputs, timeout_s=TIMEOUT['locate'])

    def rope_locate(self, group, ref, rset, cap, shift, batch='1m', figures=()):
        inputs = dict(binary=self.rope(), runtime=self.runtime(), fmd=self.input(f'ropebwt3/{ref}/index.fmd'),
                      ssa=self.input(f'ropebwt3/{ref}/index.fmd.s{shift}.ssa'),
                      lengths=self.input(f'ropebwt3/{ref}/index.fmd.len.gz'), **self.workload(ref, rset))
        return dict(id=f'{group}/cap{cap}-rope-s{shift}', group=group, figures=list(figures), tool='ropebwt3',
                    kind='locate', cap=cap, shift=shift, batch=batch, reference=ref, read_set=rset,
                    inputs=inputs, timeout_s=TIMEOUT['locate'])


def default_index(layout, ref):
    """D9: HPRC2 MS dual runs on the full + seqids index that also serves locate."""
    if layout == 'dual':
        return 'index.tms.full.seqids.orbitidx' if ref == 'human472' else 'index.tms.full.orbitidx'
    return {'ms': 'index.ms.orbitidx', 'psi': 'index.tms.psi.orbitidx', 'phi': 'index.tms.phi.orbitidx'}[layout]


def locate_index(layout):
    """D9: psi locates on the psi-locate index; phi and dual on the shared full + seqids index."""
    return 'index.tms.psi-locate.orbitidx' if layout == 'psi' else 'index.tms.full.seqids.orbitidx'


LEGACY_FILES = {'moni': ('index.thrbv.ms', 'index.plain.slp', 'index.idx'),
                'phoni': ('reference.txt.phoni', 'reference.txt.plain.slp', 'reference.txt'),
                'augphoni': ('reference.txt.aug', 'reference.txt.plain.slp', 'reference.txt')}
LEGACY_BINARY = {'moni': 'moni-ms', 'phoni': 'phoni', 'augphoni': 'aug_phoni'}
# The C++ runtime the legacy binaries link (ms2exp.build.tools), shipped for hosts that lack it.
LEGACY_RUNTIME = ('libstdc++.so.6', 'libgcc_s.so.1')


def freeze_matrix(b):
    block = 64 << 20
    # Fig 3: {mtb152, HPRC2} x {100 kb pieces, PacBio} x {min, psi, phi, dual} x {ordinary, oracle}, K = 32.
    # D13.2: the pieces are a seeded sample of 600 (60 Mb); D13.4: the PacBio sets are capped at 60,000,000 bases.
    fig3_sets = {'mtb152': ('mtb152-heldout16-p100k-s600-v1', 'mtb152-pacbio-60mb-v2'),
                 'human472': ('human472-apr008-hap1-p100k-s600-v1', 'human472-apr008-pacbio-60mb-v2')}
    for ref, sets in fig3_sets.items():
        for rset in sets:
            group = f'fig_ops-{rset}'
            cap = PHI_CAP_S if ref == 'human472' else None
            cells = {l: [b.ms2_ms(group, ref, rset, l, m, block_bytes=block, figures=['3'],
                                  cap_s=cap if (l, m) == ('phi', 'ordinary') else None)
                         for m in ('ordinary', 'oracle')] for l in LAYOUTS}
            if ref == 'human472':
                # One job per index, so each layout's cells start once its own index exists.
                for layout, layout_cells in cells.items():
                    b.job(f'fig3-{rset}-{layout}', layout_cells)
            else:
                b.job(f'fig3-{rset}', [c for l in LAYOUTS for c in cells[l]])
            counts = [b.ms2_count(group + '-counts', ref, rset, l, m, block_bytes=block, figures=['3'])
                      for l in LAYOUTS for m in ('ordinary', 'oracle')]
            b.job(f'fig3count-{rset}', counts, platforms=b.untimed_platforms, pool='untimed',
                  note='D7: counting build; deterministic, so any machine')
    # Fig 4 / A2: chr22 h100 and h10 x simulated short, PacBio, ONT, from APR008 hap1 chr22 (D13.1).
    for h in (100, 10):
        ref = f'chr22-h{h}'
        for s in CHR22_SETS:
            rset = f'{s}-60mb-v2'
            group = f'chr22-h{h}-ms-{s}'
            fig = '4' if h == 100 else 'A2'
            cells = [b.ms2_ms(group, ref, rset, l, 'ordinary', figures=[fig]) for l in LAYOUTS]
            cells += [b.rope_ms(group, ref, rset, n, batch='1m', figures=[fig]) for n in ROPE_LENGTHS]
            cells += [b.legacy(group, ref, rset, t, figures=[fig]) for t in LEGACY_FILES]
            b.job(f'fig4-{ref}-{s}', cells)
    # Fig 5 / A5: HPRC2 x 4 real 60 Mb sets.
    for cohort in HPRC2_MS_SETS:
        rset = f'{cohort}-ms-60mb-acgt-v2'
        group = f'{cohort}-ms-acgt60mb'
        for layout in LAYOUTS:
            cells = [b.ms2_ms(group, 'human472', rset, layout, m, figures=['5'],
                              cap_s=PHI_CAP_S if (layout, m) == ('phi', 'ordinary') else None)
                     for m in ('ordinary', 'oracle')]
            if cohort == 'human_short':
                cells.append(b.ms2_ms(group, 'human472', rset, layout, 'ftab', figures=['A5']))
            if layout == 'dual':
                # D9 gate: MS dual on the plain full index against dual on full + seqids, same box.
                cells.append(b.ms2_ms(group, 'human472', rset, 'dual', 'ordinary',
                                      index_name='index.tms.full.orbitidx', fig_id='dual-ordinary-plainfull',
                                      figures=['D9-check'], platforms=b.d9_platforms))
            b.job(f'fig5-{cohort}-ms2{layout}', cells)
        rope = [b.rope_ms(group, 'human472', rset, n, figures=['5']) for n in ROPE_LENGTHS]
        b.job(f'fig5-{cohort}-rope', rope)
    # Fig 6 / A3 / A4 / A6: HPRC2 locate, split by read set and index so no job runs past about 12 h.
    for cohort, size in HPRC2_LOCATE_SETS.items():
        rset = f'{cohort}-locate-{size}-v2'
        group = f'{cohort}-locate-{size}'
        figs = {1: 'A3', 10: 'A4', 100: '6'}
        # One job per layout and mode: the short-read set's phi and dual jobs would otherwise be
        # the longest of the run (about 3.6 h together).
        for layout in ('psi', 'phi', 'dual'):
            b.job(f'loc-{cohort}-ms2{layout}',
                  [b.ms2_locate(group, 'human472', rset, layout, c, figures=[figs[c]]) for c in CAPS])
            if cohort == 'human_short':
                b.job(f'loc-{cohort}-ms2{layout}-ftab',
                      [b.ms2_locate(group, 'human472', rset, layout, c, mode='ftab', figures=['A6']) for c in CAPS])
        for shift in (11, 8, 6):
            b.job(f'loc-{cohort}-rope-s{shift}',
                  [b.rope_locate(group, 'human472', rset, c, shift, figures=[figs[c]]) for c in CAPS])
    # A1: TMS_REP_LOG runs (D7), the mtb152 error series A2 and the chr22 ladder B2 (design 2 of
    # fig_rep_tikz.py). One untimed run per cell; fig_rep.py checks the logs against each other.
    b.job('a1-mtb152-errors', [b.ms2_rep('fig_rep-a2', 'mtb152', f'mtb152-A2-{tag}-v1', f'mtb152-A2-{tag}')
                               for tag in A1_ERROR_TAGS], platforms=b.untimed_platforms, pool='untimed')
    b.job('a1-chr22-ladder', [b.ms2_rep('fig_rep-b2', f'chr22-h{h}', A1_LADDER_READS.split('/')[2], f'chr22-B-h{h}',
                                        reads=A1_LADDER_READS) for h in A1_LADDER],
          platforms=b.untimed_platforms, pool='untimed')


# A1 read sets: error rates of the mtb152 series (fig_rep.py simulate, seeds 2026102801 to 10).
A1_ERROR_TAGS = ('exact', 'e0.000032', 'e0.00011', 'e0.00019', 'e0.00030', 'e0.00047', 'e0.00076', 'e0.0014',
                 'e0.0030', 'e0.0050')
# Ladder rungs of series B (D13.3); h10 and h100 are the Fig 4 / A2 indexes.
A1_LADDER = (1, 2, 5, 10, 20, 50, 100)
# The simulated APR008 chr22 PacBio reads before truth, the same reads Fig 4 uses.
A1_LADDER_READS = 'reads/chr22-sim/pacbio-60mb-v2/reads.fa'
FTAB_K = 14
# References built from a ropebwt3 FMD (one terminator per string); the A1 ladder's other rungs
# come from separator text with one terminator.
# ms2 build --layout for each index the cells name.
BUILD_LAYOUT = {'ms': 'min', 'psi': 'psi', 'phi': 'phi', 'full': 'full', 'psi-locate': 'psi-locate'}
FMD_REFS = ('human472', 'mtb152', 'chr22-h10', 'chr22-h100')


def thread_cells(b, group, rset, threads):
    """Five configurations at each thread count: ms2 min, psi, phi and dual, and ropebwt3 -l1."""
    out = {}
    for t in threads:
        binding = dict(threads=t, binding='smt' if t > 32 else 'physical', hash_output=True)
        for layout in LAYOUTS:
            cell = b.ms2_ms(group, 'human472', rset, layout, 'ordinary', block_bytes=THREAD_BLOCK,
                            figures=['threads'], fig_id=f'{layout}-ordinary-t{t}')
            out.setdefault(layout, []).append(dict(cell, **binding))
        cell = b.rope_ms(group, 'human472', rset, 1, figures=['threads'])
        out.setdefault('rope', []).append(dict(cell, id=f'{group}/rope-l1-t{t}', **binding))
    return out


def threads_matrix(b):
    """One job per read set and configuration, each running its seven thread counts on one box.
    The first job is one short cell at 64 threads, so the SMT binding is seen to work before the
    long jobs start (decision 10)."""
    smt = thread_cells(b, 'threads-smt-check', 'human_ont_modern-ms-60mb-acgt-v2', (64,))['ms']
    b.job('threads-smt-check', smt, platforms=['aws'], note='decision 10: one short 64-thread cell first')
    for cohort in THREAD_SETS:
        cells = thread_cells(b, f'threads-{cohort}', f'{cohort}-ms-1gb-acgt-v1', THREADS)
        for config, config_cells in cells.items():
            b.job(f'threads-{cohort}-{config}', config_cells, platforms=['aws'])


def threads_pilot_matrix(b):
    """Decision 10: every configuration and read set at 1 and 32 threads on the 60 Mb prefixes, one job
    per configuration. The frozen ms2 indexes are in S3 only, so the ms2 jobs run on one AWS box, after
    the 64-thread SMT cell the sweep also starts with; ropebwt3 runs on Rockfish."""
    smt = thread_cells(b, 'threads-smt-check', 'human_ont_modern-ms-60mb-acgt-v2', (64,))['ms']
    b.job('threads-smt-check', smt, platforms=['aws'], note='decision 10: one short 64-thread cell first')
    jobs = {}
    for cohort in THREAD_SETS:
        for config, config_cells in thread_cells(b, f'threads-pilot-{cohort}', f'{cohort}-ms-60mb-acgt-v2',
                                                 PILOT_THREADS).items():
            jobs.setdefault(config, []).extend(config_cells)
    for config, cells in jobs.items():
        b.job(f'threads-pilot-{config}', cells, platforms=['rockfish'] if config == 'rope' else ['aws'])


def build_matrix(b):
    """All-AWS mode: build jobs for every ms2 binary and index the cells use (D1's ms2 steps).

    Commands follow ms2's usage text at launch (639f051): `ms2 build HEADS LENS OUT --layout L`,
    where layout ms of the cells is ms2's min layout. build writes the LCP summaries itself, so
    there is no index-add-lcp-summary step. References whose BWT has one terminator per string
    (FMD_REFS) take phi from TeraLCP -ophi samples (--phi) and the summaries from the minima:
    build --minima alone walks LF over all n rows and aborts there.
    Inputs that ms2 does not produce (precursors, ropebwt3, read sets) stay pinned inputs (D1)."""
    c7 = b.ms2_commit[:7]
    tools = [f'tools/ms2-{c7}/{n}' for n in ('ms2-bench', 'ms2-stats', 'ms2-replog')] + ['tools/ms_check']
    if b.ms2_source == 'bundle':
        # langmead-lab/ms2 is private: the source comes as pinned bundles (make_source_bundles.sh).
        bundle = {m: '{in:' + f'source/ms2-{c7}/{m}.bundle' + '}' for m in ('ms2', 'orbit', 'teratools')}
        fetch = ['--ms2-repo', bundle['ms2'], '--submodule-repo', 'orbit', bundle['orbit'],
                 '--submodule-repo', 'teratools', bundle['teratools']]
    else:
        fetch = ['--ms2-repo', b.ms2_repo]
    b.build('build-tools', [dict(label='build_tools', argv=[
        'python3', '-m', 'ms2exp.build.tools', '--parts', 'ms2', *fetch, '--ms2-commit', b.ms2_commit,
        '--root', '{outroot}', '--work', '{work}/ms2-src'])], tools)
    ms = '{in:' + f'tools/ms2-{c7}/ms2-bench' + '}'

    def pre(ref, name):
        return '{in:' + f'precursors/{ref}/{name}' + '}'

    def index_steps(ref, layout, out, ftab):
        heads, lens, minima = pre(ref, 'bwt.heads'), pre(ref, 'bwt.len'), pre(ref, 'min')
        lcp = ['--phi', pre(ref, 'phi'), '--minima', minima] if ref in FMD_REFS else ['--minima', minima]
        o = '{out:' + out + '}'
        build = [ms, 'build', heads, lens, o, '--layout', BUILD_LAYOUT[layout]]
        if layout == 'ms':
            build += ['--minima', minima]
        elif layout == 'psi-locate':
            build += [*lcp, '--seq-ids', '{in:' + f'ropebwt3/{ref}/index.fmd.len.gz' + '}']
        elif layout != 'psi':
            build += lcp
        steps = [dict(label=f'build-{BUILD_LAYOUT[layout]}', argv=build)]
        outputs = [out]
        if ftab:
            steps.append(dict(label='ftab-build', argv=[ms, 'ftab-build', o, '{out:' + out + '.ftab}',
                                                        '--kmer-len', str(FTAB_K), '-p', '{threads}']))
            outputs.append(out + '.ftab')
        return steps, outputs

    # HPRC2: one box per index, so phi, full and psi-locate (about 10 h each) build side by side.
    names = {'ms': 'index.ms.orbitidx', 'psi': 'index.tms.psi.orbitidx', 'phi': 'index.tms.phi.orbitidx',
             'full': 'index.tms.full.orbitidx', 'psi-locate': 'index.tms.psi-locate.orbitidx'}
    for layout in ('ms', 'psi', 'phi', 'psi-locate'):
        steps, outputs = index_steps('human472', layout, f'ms2/human472/{names[layout]}', ftab=True)
        b.build(f'build-human472-{layout}', steps, outputs)
    plain, seqids = 'ms2/human472/index.tms.full.orbitidx', 'ms2/human472/index.tms.full.seqids.orbitidx'
    steps, outputs = index_steps('human472', 'full', plain, ftab=False)
    steps += [dict(label='index-add-seq-ids', argv=[ms, 'index-add-seq-ids', '{out:' + plain + '}',
                                                   '{in:ropebwt3/human472/index.fmd.len.gz}', '{out:' + seqids + '}']),
              dict(label='ftab-build', argv=[ms, 'ftab-build', '{out:' + seqids + '}', '{out:' + seqids + '.ftab}',
                                             '--kmer-len', str(FTAB_K), '-p', '{threads}'])]
    b.build('build-human472-full', steps, outputs + [seqids, seqids + '.ftab'])
    # Small references: every layout on one box. The A1 ladder needs psi, phi and full.
    small = {'mtb152': ('ms', 'psi', 'phi', 'full'), 'chr22-h10': ('ms', 'psi', 'phi', 'full'),
             'chr22-h100': ('ms', 'psi', 'phi', 'full')}
    for h in A1_LADDER:
        small.setdefault(f'chr22-h{h}', ('psi', 'phi', 'full'))
    for ref, layouts in small.items():
        steps, outputs = [], []
        for layout in layouts:
            st, out = index_steps(ref, layout, f'ms2/{ref}/{names[layout]}', ftab=False)
            steps += st
            outputs += out
        b.build(f'build-{ref}', steps, outputs)




def toy_matrix(b):
    """mtb152-sized cells covering every runner path, on the freeze's read sets or cuts of them; inputs come
    from toy_prep.sh."""
    block = 64 << 20
    ref = 'mtb152'
    for rset in ('mtb152-heldout16-p100k-s600-v1', 'mtb152-pacbio-60mb-v2'):
        group = f'fig_ops-{rset}'
        cells = [b.ms2_ms(group, ref, rset, l, m, block_bytes=block, figures=['3'])
                 for l in LAYOUTS for m in ('ordinary', 'oracle')]
        b.job(f'toy-fig3-{rset}', cells)
    rset = 'mtb152-pacbio-60mb-v2'
    group = f'toy-ms-{rset}'
    cells = [b.ms2_ms(group, ref, rset, l, 'ftab', figures=['A5']) for l in ('ms', 'psi')]
    cells += [b.rope_ms(group, ref, rset, n, batch='1m', figures=['4']) for n in (1, 19, 30)]
    cells.append(b.ms2_ms(group, ref, rset, 'phi', 'ordinary', fig_id='phi-ordinary-capped', cap_s=1,
                          figures=['cap-check']))
    b.job('toy-ms-mixed', cells)
    b.job('toy-count', [b.ms2_count(group + '-counts', ref, rset, l, 'ordinary', figures=['3'])
                        for l in ('ms', 'psi')], platforms=['rockfish'])
    rset = 'mtb152-locate-2mb-v2'
    group = 'toy-locate-2mb'
    cells = [b.ms2_locate(group, ref, rset, l, c, index_name='index.tms.full.seqids.orbitidx')
             for l in ('phi', 'dual') for c in (1, 10)]
    cells += [b.rope_locate(group, ref, rset, c, 8) for c in (1, 10)]
    b.job('toy-locate', cells)
    # MONI, PHONI and AugPHONI on their chr22-h10 indexes (D1), next to ms2-min on the same reads.
    ref, rset = 'chr22-h10', 'simshort-4mb-v2'
    group = f'toy-{ref}-ms-simshort'
    cells = [b.legacy(group, ref, rset, t, figures=['4']) for t in LEGACY_FILES]
    cells.append(b.ms2_ms(group, ref, rset, 'ms', 'ordinary', figures=['4']))
    b.job('toy-legacy', cells)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--matrix', choices=('freeze', 'toy', 'threads', 'threads-pilot'), required=True)
    p.add_argument('--inputs', required=True, help='inputs catalog (ms2exp.catalog hash)')
    p.add_argument('--ms2-repo', required=True, help='ms2 source repository URL')
    p.add_argument('--ms2-commit', required=True, help='full 40-hex ms2 commit')
    p.add_argument('--ms2-source', choices=('bundle', 'repo'), default='bundle',
                   help='--all-aws: build ms2 from the pinned inputs source/ms2-<c7>/{ms2,orbit,teratools}.bundle '
                        '(default; make_source_bundles.sh) or by fetching --ms2-repo')
    p.add_argument('--ms2-build', default='make -C src MS_COMMIT={commit} ms2-bench',
                   help='build command recorded for the ms2 binary (ms2exp.build.tools runs it)')
    p.add_argument('--ropebwt3-repo', default='https://github.com/lh3/ropebwt3.git')
    p.add_argument('--ropebwt3-commit', default='b408339d0fb08fe0f4da5140507a6b8b854fa12f')
    p.add_argument('--ropebwt3-binary-sha256',
                   default='cfad8475911f9298bf53177bf4a578a6d133f0464f4b9970bef9afb04697d53d',
                   help='sha256 of the pinned ropebwt3 binary of the e637ac0-era runs')
    p.add_argument('--seed', type=int, default=20261029)
    p.add_argument('--allow-pending', action='store_true')
    p.add_argument('--all-aws', action='store_true',
                   help='add build jobs and run counting, A1 and D9 cells on AWS too (no Rockfish dry run)')
    p.add_argument('--note')
    p.add_argument('-o', '--output', required=True)
    args = p.parse_args()
    if len(args.ms2_commit) != 40:
        raise SystemExit('--ms2-commit must be a full 40-hex commit')
    catalog = json.loads(Path(args.inputs).read_text())['inputs']
    b = Builder(catalog, args.allow_pending, args.ms2_commit, args.ropebwt3_binary_sha256[:8], all_aws=args.all_aws)
    b.ms2_repo = args.ms2_repo
    b.ms2_source = args.ms2_source
    if args.all_aws:
        build_matrix(b)
    matrices = {'freeze': freeze_matrix, 'toy': toy_matrix, 'threads': threads_matrix,
                'threads-pilot': threads_pilot_matrix}
    matrices[args.matrix](b)
    if args.all_aws:
        unbuilt = sorted(n for n in b.inputs if n.startswith(('ms2/', 'tools/ms2-')) and n not in b.produced)
        if unbuilt:
            raise SystemExit('No build job produces:\n  ' + '\n  '.join(unbuilt))
    if b.missing:
        raise SystemExit('Inputs missing from the catalog (or pass --allow-pending):\n  ' + '\n  '.join(b.missing))
    rope_binary = b.inputs.get(f'tools/ropebwt3-{args.ropebwt3_binary_sha256[:8]}/ropebwt3', {})
    if rope_binary.get('sha256') not in ('pending', args.ropebwt3_binary_sha256):
        raise SystemExit('The cataloged ropebwt3 binary differs from --ropebwt3-binary-sha256')
    manifest = dict(
        format=MANIFEST_FORMAT, matrix=args.matrix, created=now(), generator='ms2exp/manifest.py',
        note=args.note, seed=args.seed, protocol=THREADS_PROTOCOL if args.matrix.startswith('threads') else PROTOCOL,
        ms2=dict(repo=args.ms2_repo, commit=args.ms2_commit, build=args.ms2_build.format(commit=args.ms2_commit),
                 **(dict(source=args.ms2_source) if args.all_aws else {})),
        ropebwt3=dict(repo=args.ropebwt3_repo, commit=args.ropebwt3_commit,
                      binary_sha256=args.ropebwt3_binary_sha256),
        inputs=dict(sorted(b.inputs.items())), jobs=b.jobs)
    manifest['id'] = f'{args.matrix}-{args.ms2_commit[:7]}-{json_sha256(manifest)[:8]}'
    Path(args.output).write_text(json.dumps(manifest, indent=1) + '\n')
    cells = sum(len(j['cells']) for j in b.jobs)
    pending = sum(1 for v in b.inputs.values() if v['sha256'] == 'pending')
    print(f'{manifest["id"]}: {len(b.jobs)} jobs, {cells} cells, {len(b.inputs)} inputs ({pending} pending)')


if __name__ == '__main__':
    main()
