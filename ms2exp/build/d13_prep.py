#!/usr/bin/env python3
"""Read sets and truth for design decision D13 (BENCH_FREEZE_DESIGN.md, dataset consolidation).

Every set lands in OUT/<name>/ (--out; default DATA/prepared, DATA being the data directory,
MS2EXP_DATA) with reads.fa,
offsets.u64 (read boundaries, n + 1 little-endian u64), truth.u32 (per-position MS, when the
set has truth), one.fa (its first read, for the legacy one-read run) where the source has one,
and manifest.json recording the source, rule, counts and file checksums.

Steps (run on Rockfish; the heavy ones in Slurm):
  apr008-chr22 --scaffold NAME
      D13.1: the APR008 hap1 scaffold that is chr22 (found by aligning every scaffold to CHM13
      chr22, apr008-chr22/per-scaffold.tsv), uppercased and split at non-ACGT into ACGT-only
      records: OUT/haplotype-APR008-hap1-chr22-v1/acgt.fa.
  simulate --set short|pacbio|ont
      D13.1 and D13.4: chr22_sim.py's recipe (simulator, model, settings, seed) on that source,
      then the whole-read prefix of at most 60,000,000 bases, then ropebwt3 -t1 -l1 truth on the
      chr22-h10 and chr22-h100 FMDs (chr22_sim.truth). Writes chr22-sim-SET-60mb-v2 (reads only,
      the A1 ladder reads) and chr22-hN-SET-60mb-v2 (with truth).
  cap SRC NAME --limit BASES
      D13.4 and D13.5: the longest whole-read prefix of SRC with at most BASES bases; truth,
      offsets and one.fa are cut to match. SRC's reads must already be ACGT-only whole reads.
  sample SRC NAME --count 600 --seed S
      D13.2: a seeded random sample of SRC's records (kept in source order), truth sliced.
  short-60mb
      D13.4 for the real short reads: the historical core pool (1.6 M records, 199.5 Mb, with
      validated truth) filtered first to reads that were whole and ACGT-only, then the first
      reads in core order up to 60,000,000 bases; truth sliced from the core truth.
  skin-60mb [--threads N]
      D13.4 for the skin reads: ERR14066769 pairs in deposited order, R1 then R2, uppercased;
      reads with any non-ACGT base dropped; the first reads up to 60,000,000 bases; truth by
      ropebwt3 -l1 on the human472 FMD, converted as the historical skin truth was.
  superset --cohort human_short|human_ont_modern NAME --prefix SET [--limit BASES]
      Thread-scaling read sets (THREAD_SCALING_PLAN.md decision 4): the pool the 60 Mb set was
      drawn from (real-hprc2-1d580b6 selected.fa) under the D13.4 filter, its core records first
      in core order and then its other records in pool order, stopping at the last read that
      keeps the total at or below BASES (default 1,000,000,000). The 60 Mb set SET is a prefix of
      the core records in core order, so it must be a byte-exact prefix of the result (reads.fa
      and offsets.u64), which this step checks. Reads only; the truth step adds truth.u32.
  truth NAME [--threads N] [--work DIR]
      ropebwt3 cfad8475 -l1 truth on the human472 FMD for a prepared set without truth,
      converted and bounds-checked as the skin truth was; when the set names a prefix set, the
      new truth over that set's positions must equal the prefix set's truth byte for byte.

A set is never overwritten: a step whose destination exists stops.

Steps that start from lab-internal pools (apr008-chr22's raw assembly, short-60mb, skin-60mb's
FASTQ files and superset) read them under the lab root, MS2EXP_LAB_ROOT. ropebwt3 and the human472
FMD come from the data directory by logical name (tools/ropebwt3-cfad8475/ropebwt3,
ropebwt3/human472/index.fmd), or from --ropebwt3 and --fmd, and the binary is checked against its pin.
"""
import argparse
import array
import gzip
import json
import os
from pathlib import Path
import random
import re
import shutil
import sys
import time

from ms2exp import paths
from ms2exp.run import common
from ms2exp.sim import chr22_sim
from ms2exp.sim.chr22_sim import sha256
from ms2exp.validate.ms import check

ROOT = paths.LAB_ROOT

P = ROOT/'results/campaign-20261012'
OUT = None
LIMIT = 60_000_000
FILTER = 'acgt-only-whole-reads'
ACGT = re.compile('^[ACGT]+$')
APR008_RAW = P/'inputs/heldout/apr-APR008-hap1/raw/GCA_050494745.1_apr008.1_v1_genomic.fna.gz'
SHORT_CORE = ROOT/'results/real-hprc2-1d580b6/inputs/human_short'
SHORT_TRUTH = ROOT/'results/real-hprc2-1d580b6/truth/human_short/core.u32'
SKIN = ROOT/'reads/guccione-PRJEB83637/skin/ERR14066769'
POOLS = {c: ROOT/'results/real-hprc2-1d580b6/inputs'/c for c in ('human_short', 'human_ont_modern')}
SUPERSET_LIMIT = 1_000_000_000
ROPE_SHA256 = 'cfad8475911f9298bf53177bf4a578a6d133f0464f4b9970bef9afb04697d53d'


def rope_and_fmd(args):
    """The pinned ropebwt3 and the human472 FMD: --ropebwt3 and --fmd, else the data directory's."""
    data = None if args.ropebwt3 and args.fmd else paths.data_dir()
    rope = Path(args.ropebwt3 or data/'tools/ropebwt3-cfad8475/ropebwt3')
    fmd = Path(args.fmd or data/'ropebwt3/human472/index.fmd')
    if sha256(rope) != ROPE_SHA256:
        raise ValueError('ropebwt3 differs from the frozen comparator')
    return rope, fmd


def read_u64(path):
    a = array.array('Q')
    a.frombytes(Path(path).read_bytes())
    if sys.byteorder != 'little':
        a.byteswap()
    return a


def write_u64(path, values):
    a = array.array('Q', values)
    if sys.byteorder != 'little':
        a.byteswap()
    Path(path).write_bytes(a.tobytes())


def records(path):
    """(name, sequence) of a FASTA whose records are one header line and one sequence line."""
    with Path(path).open() as f:
        while True:
            head = f.readline()
            if not head:
                return
            yield head[1:].rstrip('\n').split()[0], f.readline().rstrip('\n')


def finish(dest, manifest):
    files = sorted(p.name for p in dest.iterdir() if p.is_file() and p.name != 'manifest.json')
    manifest['files'] = {n: dict(bytes=(dest/n).stat().st_size, sha256=sha256(dest/n)) for n in files}
    manifest.setdefault('created', time.time())
    manifest.setdefault('job_id', os.environ.get('SLURM_JOB_ID'))
    (dest/'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
    print(json.dumps({k: manifest[k] for k in ('name', 'reads', 'bases') if k in manifest}), flush=True)
    return manifest


def new_dir(name):
    dest = OUT/name
    dest.mkdir(parents=True, exist_ok=False)
    return dest


def source_record(src):
    src = Path(src)
    m = src/'manifest.json'
    return dict(path=str(src), manifest_sha256=sha256(m) if m.exists() else None)


def write_subset(src, dest, chosen, offsets):
    """Write the records of src with the given indices (ascending), with offsets and sliced truth."""
    chosen = list(chosen)
    keep = set(chosen)
    new = [0]
    with (dest/'reads.fa').open('x') as out:
        for i, (_, seq) in enumerate(records(src/'reads.fa')):
            if i in keep:
                if len(seq) != offsets[i + 1] - offsets[i]:
                    raise ValueError(f'Read {i} length differs from offsets')
                out.write(f'>q{len(new) - 1}\n{seq}\n')
                new.append(new[-1] + len(seq))
            if i >= chosen[-1]:
                break
    write_u64(dest/'offsets.u64', new)
    if (src/'truth.u32').exists():
        with (src/'truth.u32').open('rb') as t, (dest/'truth.u32').open('xb') as out:
            for i in chosen:
                t.seek(4 * offsets[i])
                out.write(t.read(4 * (offsets[i + 1] - offsets[i])))
    return new


def cmd_cap(args):
    src = Path(args.src)
    offsets = read_u64(src/'offsets.u64')
    k = max(i for i in range(len(offsets)) if offsets[i] <= args.limit)
    if k == 0:
        raise ValueError('No whole read fits the limit')
    dest = new_dir(args.name)
    new = write_subset(src, dest, range(k), offsets)
    if (src/'one.fa').exists():
        shutil.copyfile(src/'one.fa', dest/'one.fa')
    for table in ('origins.tsv', 'map.tsv'):
        if (src/table).exists():
            with (src/table).open() as f, (dest/table).open('x') as out:
                lines = f.readlines()
                header = 1 if lines and not lines[0].startswith('q') else 0
                out.writelines(lines[:header + k])
    old = json.loads((src/'manifest.json').read_text()) if (src/'manifest.json').exists() else {}
    manifest = dict(old, name=args.name, status=old.get('status', 'prepared'), reads=k, bases=new[-1],
                    read_filter=old.get('read_filter', FILTER), limit_bases=args.limit,
                    selection=f'longest whole-read prefix of the source with at most {args.limit} bases (D13)',
                    source=source_record(src), source_reads=len(offsets) - 1, source_bases=offsets[-1])
    for stale in ('files', 'target_bases', 'lengths', 'strands', 'read_lengths'):
        manifest.pop(stale, None)
    finish(dest, manifest)


def cmd_sample(args):
    src = Path(args.src)
    offsets = read_u64(src/'offsets.u64')
    n = len(offsets) - 1
    chosen = sorted(random.Random(args.seed).sample(range(n), args.count))
    dest = new_dir(args.name)
    new = write_subset(src, dest, chosen, offsets)
    with (dest/'map.tsv').open('x') as out:
        out.write('name\tsource_record\n')
        for j, i in enumerate(chosen):
            out.write(f'q{j}\tq{i}\n')
    old = json.loads((src/'manifest.json').read_text())
    manifest = dict(name=args.name, status='prepared', reference=old.get('reference'), cohort=old.get('cohort'),
                    reads=len(chosen), bases=new[-1], read_filter=old.get('read_filter'),
                    selection=(f'seeded random sample (Python random.Random({args.seed}).sample) of {args.count} of '
                               f'the {n} source records, kept in source order and renamed q0, q1, ...; map.tsv '
                               'names each source record (D13)'),
                    seed=args.seed, source=source_record(src), source_reads=n, source_bases=offsets[-1],
                    truth='sliced from the source truth.u32' if (src/'truth.u32').exists() else None)
    finish(dest, manifest)


def cmd_apr008_chr22(args):
    dest = new_dir('haplotype-APR008-hap1-chr22-v1')
    seq = None
    with gzip.open(APR008_RAW, 'rt') as f:
        name, parts = None, []
        for line in f:
            if line.startswith('>'):
                if name == args.scaffold:
                    break
                name, parts = line[1:].split()[0], []
            elif name == args.scaffold:
                parts.append(line.strip())
        if name == args.scaffold:
            seq = ''.join(parts).upper()
    if not seq:
        raise ValueError('Scaffold not found: ' + args.scaffold)
    segments = [(m.start(), m.group()) for m in re.finditer('[ACGT]+', seq)]
    with (dest/'acgt.fa').open('x') as out:
        for i, (start, s) in enumerate(segments):
            out.write(f'>APR008_hap1_chr22_s{i}_{start}\n{s}\n')
    manifest = dict(name=dest.name, format='chr22-sim-v1', status='prepared', sample='APR008 hap1',
                    accession='GCA_050494745.1', scaffold=args.scaffold, scaffold_length=len(seq),
                    acgt_bases=sum(len(s) for _, s in segments), acgt_segments=len(segments),
                    non_acgt_bases=len(seq) - sum(len(s) for _, s in segments), raw=str(APR008_RAW),
                    raw_sha256=sha256(APR008_RAW), identification=args.evidence,
                    selection='the scaffold uppercased and split at every non-ACGT character into ACGT-only records')
    finish(dest, manifest)


def cmd_simulate(args):
    hap = OUT/'haplotype-APR008-hap1-chr22-v1'
    note = f'APR008 hap1 chr22 ({json.loads((hap/"manifest.json").read_text())["scaffold"]}), D13'
    source = chr22_sim.load_source(hap/'acgt.fa', note=note)
    raw = OUT/f'chr22-sim-{args.set}-raw-v2'
    if not raw.exists():
        chr22_sim.simulate(args.set, threads=args.threads, destination=raw, source=source,
                           cohort=chr22_sim.SETS[args.set]['cohort'])
    reads = OUT/f'chr22-sim-{args.set}-60mb-v2'
    if not reads.exists():
        cmd_cap(argparse.Namespace(src=raw, name=reads.name, limit=LIMIT))
        # chr22_sim.truth reads this set's manifest as prepared simulated reads.
        m = json.loads((reads/'manifest.json').read_text())
        m.update(status='reads-prepared', set=args.set, read_filter=chr22_sim.FILTER)
        finish(reads, m)
    for reference in chr22_sim.REFERENCES:
        dest = OUT/f'{reference}-sim{args.set}-60mb-v2'
        if not dest.exists():
            chr22_sim.truth(args.set, reference, source=reads, destination=dest,
                            work=OUT.parent/'truth-jobs'/f'{reference}-{args.set}-{os.environ.get("SLURM_JOB_ID", "local")}')


def cmd_short(args):
    names = {}
    with (SHORT_CORE/'original-read-names.tsv').open() as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            names[fields[0]] = int(fields[-1])
    offsets = read_u64(SHORT_CORE/'core.offsets.u64')
    rows = []
    with (SHORT_CORE/'core.map.tsv').open() as f:
        for i, line in enumerate(f):
            q, key, start, end = line.split('\t')[:4]
            if q != f'q{i}':
                raise ValueError('Core map order differs')
            rows.append((key, int(start), int(end)))
    if len(rows) != len(offsets) - 1:
        raise ValueError('Core map and offsets differ')
    count = {}
    for key, _, _ in rows:
        count[key] = count.get(key, 0) + 1
    chosen, bases = [], 0
    for i, (key, start, end) in enumerate(rows):
        if count[key] != 1 or start != 0 or end != names[key]:
            continue
        n = offsets[i + 1] - offsets[i]
        if bases + n > LIMIT:
            break
        chosen.append(i)
        bases += n
    dest = new_dir('human_short-ms-60mb-acgt-v2')
    src = SHORT_CORE
    # write_subset reads src/reads.fa and src/truth.u32; the core keeps them under other names.
    stage = OUT.parent/'stage-short'
    stage.mkdir(parents=True, exist_ok=True)
    for a, b in (('core.fa', 'reads.fa'),):
        if not (stage/b).exists():
            os.symlink(src/a, stage/b)
    if not (stage/'truth.u32').exists():
        os.symlink(SHORT_TRUTH, stage/'truth.u32')
    new = write_subset(stage, dest, chosen, offsets)
    (dest/'one.fa').write_text(f">q0\n{next(records(dest/'reads.fa'))[1]}\n")
    with (dest/'map.tsv').open('x') as out:
        out.write('name\tcore_record\toriginal_read\n')
        for j, i in enumerate(chosen):
            out.write(f'q{j}\tq{i}\t{rows[i][0]}\n')
    manifest = dict(name=dest.name, status='prepared', reference='human472', cohort='human_short', reads=len(chosen),
                    bases=new[-1], one_bases=new[1], read_filter=FILTER, limit_bases=LIMIT,
                    selection=('core records whose original read became exactly one record covering it whole (so it '
                               'was entirely A, C, G and T after uppercasing), in core order, stopping at the last read '
                               f'that keeps the total at or below {LIMIT} bases (D13.4)'),
                    source=dict(core=str(SHORT_CORE/'core.fa'), core_sha256=sha256(SHORT_CORE/'core.fa'),
                                truth=str(SHORT_TRUTH), truth_sha256=sha256(SHORT_TRUTH)),
                    truth='sliced from the core truth')
    finish(dest, manifest)


def fastq(path):
    with gzip.open(path, 'rt') as f:
        while True:
            head = f.readline()
            if not head:
                return
            seq = f.readline().rstrip('\n')
            f.readline()
            f.readline()
            yield head[1:].split()[0], seq


def cmd_skin(args):
    dest = new_dir('guccione_skin-ms-60mb-acgt-v2')
    m1, m2 = SKIN/'ERR14066769_1.fastq.gz', SKIN/'ERR14066769_2.fastq.gz'
    offsets, bases, dropped = [0], 0, 0
    done = False
    with (dest/'reads.fa').open('x') as out, (dest/'map.tsv').open('x') as table:
        table.write('name\tpair_ordinal\tmate\tdeposited_name\n')
        for ordinal, ((n1, s1), (n2, s2)) in enumerate(zip(fastq(m1), fastq(m2))):
            for mate, name, seq in ((1, n1, s1), (2, n2, s2)):
                seq = seq.upper()
                if not ACGT.match(seq):
                    dropped += 1
                    continue
                if bases + len(seq) > LIMIT:
                    done = True
                    break
                q = len(offsets) - 1
                out.write(f'>q{q}\n{seq}\n')
                table.write(f'q{q}\t{ordinal}\t{mate}\t{name}\n')
                bases += len(seq)
                offsets.append(bases)
            if done:
                break
    write_u64(dest/'offsets.u64', offsets)
    (dest/'one.fa').write_text(f">q0\n{next(records(dest/'reads.fa'))[1]}\n")
    rope, fmd = rope_and_fmd(args)
    work = OUT.parent/'truth-jobs'/f"skin-{os.environ.get('SLURM_JOB_ID', 'local')}"
    work.mkdir(parents=True, exist_ok=False)
    row = common.execute([rope, 'mem', '-t', str(args.threads), '-l1', fmd, dest/'reads.fa'], work, 'truth', 6 * 3600)
    partial = work/'truth.u32'
    row['check'] = check('ropebwt3', work, row, partial, offsets, write=True)
    chr22_sim.bounded(partial, offsets)
    Path(row['stdout']).unlink()
    row['raw_output_retained'] = False
    shutil.move(str(partial), dest/'truth.u32')
    (dest/'truth.json').write_text(json.dumps(dict(passed=True, command=row, fmd=str(fmd),
                                                   ropebwt3=dict(path=str(rope), sha256=ROPE_SHA256)), indent=1, default=str))
    manifest = dict(name=dest.name, status='prepared', reference='human472', cohort='guccione_skin',
                    reads=len(offsets) - 1, bases=bases, one_bases=offsets[1], read_filter=FILTER, limit_bases=LIMIT,
                    run='ERR14066769', dropped_non_acgt_reads=dropped,
                    selection=('ERR14066769 pairs in deposited order, R1 then R2; each read uppercased and dropped if it '
                               'holds any non-ACGT base; the rest taken whole in that order, stopping at the last read '
                               f'that keeps the total at or below {LIMIT} bases (D13.4); map.tsv names each read'),
                    sources=[dict(path=str(p), sha256=sha256(p)) for p in (m1, m2)],
                    truth='ropebwt3 cfad8475 mem -l1 on the human472 FMD, converted to per-position MS by '
                          'real-hprc2 benchmark.check, every value bounds-checked')
    finish(dest, manifest)


def same_prefix(path, prefix):
    """True when the file at prefix is a byte-exact prefix of the file at path."""
    with Path(prefix).open('rb') as a, Path(path).open('rb') as b:
        while chunk := a.read(1 << 24):
            if b.read(len(chunk)) != chunk:
                return False
    return True


def cmd_superset(args):
    pool = POOLS[args.cohort]
    prefix = OUT/args.prefix if not Path(args.prefix).is_absolute() else Path(args.prefix)
    names = {}
    with (pool/'original-read-names.tsv').open() as f:
        for line in f:
            key, _, original, length = line.rstrip('\n').split('\t')
            names[key] = (original, int(length))
    # Core records in core order; a read split at non-ACGT characters appears in several records.
    core = {}
    with (pool/'core.map.tsv').open() as f:
        for line in f:
            core.setdefault(line.split('\t')[1], len(core))
    order = [int(k[1:]) for k in core]
    if order != sorted(order):
        raise ValueError('Core order differs from pool order')
    dest = new_dir(args.name)
    offsets, counts = [0], dict(core=0, rest=0)
    dropped = dict(core=0, rest=0)
    full = False
    with (dest/'reads.fa').open('x') as out, (dest/'map.tsv').open('x') as table:
        table.write('name\tpool_record\toriginal_read\tpart\n')
        for part in ('core', 'rest'):
            for key, seq in records(pool/'selected.fa'):
                if (key in core) != (part == 'core'):
                    continue
                if len(seq) != names[key][1]:
                    raise ValueError(f'{key} length differs from original-read-names.tsv')
                seq = seq.upper()
                if not ACGT.match(seq):
                    dropped[part] += 1
                    continue
                if offsets[-1] + len(seq) > args.limit:
                    full = True
                    break
                q = len(offsets) - 1
                out.write(f'>q{q}\n{seq}\n')
                table.write(f'q{q}\t{key}\t{names[key][0]}\t{part}\n')
                offsets.append(offsets[-1] + len(seq))
                counts[part] += 1
            if full:
                break
    write_u64(dest/'offsets.u64', offsets)
    (dest/'one.fa').write_text(f">q0\n{next(records(dest/'reads.fa'))[1]}\n")
    for f in ('reads.fa', 'offsets.u64', 'one.fa'):
        if not same_prefix(dest/f, prefix/f):
            raise ValueError(f'{prefix.name}/{f} is not a prefix of {dest.name}/{f}')
    prefix_manifest = json.loads((prefix/'manifest.json').read_text())
    manifest = dict(name=dest.name, status='reads-prepared', reference='human472', cohort=args.cohort,
                    reads=len(offsets) - 1, bases=offsets[-1], one_bases=offsets[1], read_filter=FILTER,
                    limit_bases=args.limit, reached_limit=full, core_reads=counts['core'], rest_reads=counts['rest'],
                    dropped_non_acgt_reads=dropped,
                    selection=('the reads of the pool selected.fa whose original read is entirely A, C, G and T after '
                               'uppercasing (D13.4), the core records first in core order and then the other records '
                               'in pool order, taken whole and stopping at the last read that keeps the total at or '
                               f'below {args.limit} bases (THREAD_SCALING_PLAN.md decision 4); map.tsv names each read'),
                    pool=dict(path=str(pool/'selected.fa'), sha256=sha256(pool/'selected.fa'),
                              core_map=str(pool/'core.map.tsv'), core_map_sha256=sha256(pool/'core.map.tsv'),
                              records=len(names)),
                    prefix_of=dict(name=prefix.name, path=str(prefix), reads=prefix_manifest['reads'],
                                   bases=prefix_manifest['bases'], files_equal_prefix=['reads.fa', 'offsets.u64', 'one.fa']))
    finish(dest, manifest)


def cmd_truth(args):
    dest = OUT/args.name
    manifest = json.loads((dest/'manifest.json').read_text())
    if (dest/'truth.u32').exists():
        raise ValueError(f'{dest.name} already has truth')
    offsets = read_u64(dest/'offsets.u64')
    rope, fmd = rope_and_fmd(args)
    work = Path(args.work or OUT.parent/'truth-jobs')/f"{dest.name}-{os.environ.get('SLURM_JOB_ID', 'local')}"
    work.mkdir(parents=True, exist_ok=False)
    row = common.execute([rope, 'mem', '-t', str(args.threads), '-l1', fmd, dest/'reads.fa'], work, 'truth', 12 * 3600)
    partial = work/'truth.u32'
    row['check'] = check('ropebwt3', work, row, partial, offsets, write=True)
    chr22_sim.bounded(partial, offsets)
    row['stdout_bytes'] = Path(row['stdout']).stat().st_size
    Path(row['stdout']).unlink()
    row['raw_output_retained'] = False
    record = dict(passed=True, command=row, fmd=str(fmd), fmd_sha256=sha256(fmd),
                  ropebwt3=dict(path=str(rope), sha256=ROPE_SHA256))
    prefix = manifest.get('prefix_of')
    if prefix:
        theirs = Path(prefix['path'])/'truth.u32'
        if theirs.stat().st_size != 4 * prefix['bases'] or not same_prefix(partial, theirs):
            raise ValueError(f'Truth over the first {prefix["bases"]} positions differs from {prefix["name"]}')
        record['prefix_check'] = dict(name=prefix['name'], truth_sha256=sha256(theirs), bytes=4 * prefix['bases'],
                                      equal=True)
    shutil.move(str(partial), dest/'truth.u32')
    (dest/'truth.json').write_text(json.dumps(record, indent=1, default=str))
    manifest.update(status='prepared',
                    truth='ropebwt3 cfad8475 mem -l1 on the human472 FMD, converted to per-position MS by '
                          'real-hprc2 benchmark.check, every value bounds-checked' +
                          (f'; equal to {prefix["name"]} truth over its positions' if prefix else ''))
    manifest.pop('files', None)
    finish(dest, manifest)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--out', help='directory of the prepared sets (default DATA/prepared)')
    sub = p.add_subparsers(dest='step', required=True)
    a = sub.add_parser('apr008-chr22')
    a.add_argument('--scaffold', required=True)
    a.add_argument('--evidence', required=True, help='how the scaffold was identified as chr22')
    s = sub.add_parser('simulate')
    s.add_argument('--set', required=True, choices=sorted(chr22_sim.SETS))
    s.add_argument('--threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 8)))
    c = sub.add_parser('cap')
    c.add_argument('src')
    c.add_argument('name')
    c.add_argument('--limit', type=int, default=LIMIT)
    m = sub.add_parser('sample')
    m.add_argument('src')
    m.add_argument('name')
    m.add_argument('--count', type=int, default=600)
    m.add_argument('--seed', type=int, required=True)
    sub.add_parser('short-60mb')
    k = sub.add_parser('skin-60mb')
    k.add_argument('--threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 16)))
    k.add_argument('--ropebwt3', help='the pinned ropebwt3 binary (default: the data directory\'s)')
    k.add_argument('--fmd', help='the human472 FMD (default: the data directory\'s)')
    u = sub.add_parser('superset')
    u.add_argument('name')
    u.add_argument('--cohort', required=True, choices=sorted(POOLS))
    u.add_argument('--prefix', required=True, help='the 60 Mb set that must be a prefix (a name under --out, or a path)')
    u.add_argument('--limit', type=int, default=SUPERSET_LIMIT)
    t = sub.add_parser('truth')
    t.add_argument('name')
    t.add_argument('--threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 16)))
    t.add_argument('--work', help='directory for the ropebwt3 output (default OUT/../truth-jobs)')
    t.add_argument('--ropebwt3', help='the pinned ropebwt3 binary (default: the data directory\'s)')
    t.add_argument('--fmd', help='the human472 FMD (default: the data directory\'s)')
    args = p.parse_args()
    global OUT
    OUT = Path(args.out) if args.out else paths.data_dir()/'prepared'
    {'apr008-chr22': cmd_apr008_chr22, 'simulate': cmd_simulate, 'cap': cmd_cap, 'sample': cmd_sample,
     'short-60mb': cmd_short, 'skin-60mb': cmd_skin, 'superset': cmd_superset, 'truth': cmd_truth}[args.step](args)


if __name__ == '__main__':
    main()
