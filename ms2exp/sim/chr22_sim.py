"""Simulated read sets and their MS truth, as used for the chr22 and D13 read sets and Fig rep.

simulate() makes about TARGET bases of whole reads, both strands, with fixed seeds, from ACGT-only
source records (the held-out NA19240 paternal chr22 haplotype, or another FASTA given as source):
    short: ART 2.5.8 (art_illumina -ss HS25, HiSeq 2500 built-in profile, single-end 125 bases).
    pacbio: pbsim3 3.0.5 multi-pass simulation, the documented HiFi route: template strategy,
      QSHMM-RSII.model, --pass-num 10, then PacBio ccs 6.4.0 with its defaults (at least 3
      passes, predicted accuracy at least 0.99). Template lengths are drawn with replacement
      from the read lengths of the real PacBio HiFi workload and each template's strand is drawn
      with probability 1/2.
    ont: pbsim3 3.0.5 with its standard ONT error model, ERRHMM-ONT-HQ.model (--method errhmm),
      template strategy, template lengths drawn from the real HG002 Q28 Nanopore workload.
      --accuracy-mean 0.96 is the model's ceiling; higher requests give the same reads. The
      achieved accuracy and error rates pbsim3 reports are recorded.
Reads are uppercased, any read with a non-ACGT base is dropped, the rest are shuffled with the
set's seed and whole reads are taken until at least TARGET bases. They are renumbered q0, q1, ...;
origins.tsv records each read's origin and manifest.json the simulator, binaries, commands, seeds,
counts and length quantiles.

truth() runs ropebwt3 cfad8475 `mem -t1 -l1` on the reference's FMD, converts its SMEMs to
per-position MS with ms2exp.validate.ms.check, bounds-checks them (1 <= ms[j] <= length - j) and
writes them with copies of the reads as a workload.
"""
import array
import bisect
import collections
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import re
import shutil
import socket
import subprocess
import sys
import time

from ms2exp.run import common
from ms2exp.run.common import atomic_save
from ms2exp.validate.ms import check

# The ropebwt3 binary every comparison pins.
ROPE_SHA256 = 'cfad8475911f9298bf53177bf4a578a6d133f0464f4b9970bef9afb04697d53d'

FORMAT = 'chr22-sim-v1'
ROOT = Path('/scratch16/blangme2/ms-bakery')
P = ROOT/'results/campaign-20261012'
SIM = P/'inputs/chr22_sim'
FILTER = 'acgt-only-whole-reads'
TARGET = 60_000_000
REFERENCES = ('chr22-h10', 'chr22-h100')
HAPLOTYPE = SIM/'haplotype-NA19240-pat-chr22-v1'
TOOLS = ROOT/'tools'
SIM_ENV = TOOLS/'sim-env'
HIFI_ENV = TOOLS/'hifi-env/.pixi/envs/default'
BINARIES = {
    'art_illumina': (SIM_ENV/'bin/art_illumina', 'a3d0c556d5c938079d616a939325fbed37eb4b1c26413304dde550fe79ad50b4',
                     'ART_Illumina 2.5.8 (bioconda art 2016.06.05 h0704011_13)'),
    'pbsim': (SIM_ENV/'bin/pbsim', 'c58e20aaad3a39dc27a02f5840dc6eca16531e796ff79cb829c6499cedaecc30',
              'pbsim3 3.0.5 (bioconda pbsim3 3.0.5 h9948957_2)'),
    'ccs': (HIFI_ENV/'bin/ccs', '242bd8b27185ffb6d6528556197d7c4bc04bf9dc3d04149a025f587fd138b9d8',
            'ccs 6.4.0 (bioconda pbccs 6.4.0 h9ee0642_0)'),
    'samtools': (HIFI_ENV/'bin/samtools', 'd82441443bc73aad354d54aed51d6ca98a67965c992911d96c8726b1da6a140d',
                 'samtools 1.24 (bioconda samtools 1.24 h9dcdb79_1), called by pbsim3 to write BAM'),
}
MODEL = {'qshmm-rsii': SIM_ENV/'data/QSHMM-RSII.model', 'errhmm-ont-hq': SIM_ENV/'data/ERRHMM-ONT-HQ.model'}
REAL = {'short': P/'inputs/human_short/ms-60mb-acgt-v1',
        'pacbio': P/'inputs/human_pacbio/explore-acgt-60mb-v1',
        'ont': P/'inputs/human_ont_modern/explore-acgt-60mb-v1'}
SETS = {
    'short': dict(cohort='chr22_sim_short', seed=2026102501, oversample=1.10,
                  platform='Illumina short reads, ART HiSeq 2500 profile, single-end 125 bases'),
    'pacbio': dict(cohort='chr22_sim_pacbio', seed=2026102502, oversample=1.15,
                   platform='PacBio HiFi, pbsim3 multi-pass (QSHMM-RSII, 10 passes) and ccs'),
    'ont': dict(cohort='chr22_sim_ont', seed=2026102503, oversample=1.15,
                platform='Oxford Nanopore, pbsim3 ERRHMM-ONT-HQ error model, template lengths from the HG002 Q28 workload'),
}
SHORT_LENGTH = 125


# ERRHMM-ONT-HQ's ceiling: higher requests give the same reads.
ONT_ACCURACY_MEAN = 0.96
HIFI_PASSES = 10


# Simulator timeouts in seconds; ccs dominates (about 20 CPU seconds per template Mbase in a dev test).
TIMEOUT = 6 * 3600
COMPLEMENT = str.maketrans('ACGT', 'TGCA')
ACGT = frozenset('ACGT')


def reads_dir(name):
    return SIM/f'{name}-60mb-v1'


def workload_dir(reference, name):
    return SIM/f'{reference}-{name}-60mb-v1'


def verify(path, expected, what):
    if sha256(path) != expected:
        raise ValueError(f'{what} SHA256 differs: {path}')


def save_new(path, value):
    """Write JSON to a path that must not exist yet."""
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    atomic_save(path, value)


def fasta(path):
    """(name, sequence) records of a FASTA file, sequences joined across lines."""
    name, chunks = None, []
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as stream:
        for line in stream:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, ''.join(chunks)


def fastq(path):
    """(name, sequence, quality) records of a four-line FASTQ, plain or gzipped."""
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as stream:
        while True:
            header = stream.readline()
            if not header:
                return
            seq, plus, qual = stream.readline().rstrip('\n'), stream.readline(), stream.readline().rstrip('\n')
            if not header.startswith('@') or not plus.startswith('+') or len(seq) != len(qual):
                raise ValueError('Malformed FASTQ record in ' + str(path))
            yield header[1:].split()[0], seq, qual


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_u64(path):
    values = array.array('Q')
    values.frombytes(Path(path).read_bytes())
    if sys.byteorder != 'little':
        values.byteswap()
    return values


def quantiles(lengths):
    values = sorted(lengths)
    n = len(values)
    if not n:
        return {}
    pick = lambda f: values[min(n - 1, int(f * n))]
    total, running, n50 = sum(values), 0, None
    for v in reversed(values):
        running += v
        if running * 2 >= total:
            n50 = v
            break
    return dict(reads=n, bases=total, mean=total / n, min=values[0], q01=pick(.01), q10=pick(.1), q25=pick(.25),
                median=pick(.5), q75=pick(.75), q90=pick(.9), q99=pick(.99), max=values[-1], n50=n50)


def binary(name):
    path, expected, version = BINARIES[name]
    verify(path, expected, name)
    return dict(path=str(path), sha256=expected, version=version)


def load_haplotype(path=HAPLOTYPE):
    path = Path(path)
    manifest = json.loads((path/'manifest.json').read_text())
    for name, record in manifest['files'].items():
        verify(path/name, record['sha256'], name)
    segments = list(fasta(path/'acgt.fa'))
    if sum(len(s) for _, s in segments) != manifest['acgt_bases']:
        raise ValueError('ACGT segments differ from the haplotype manifest')
    return manifest, segments


def load_source(path, pair_reverse_complements=False, note=None):
    """ACGT-only segments of a source FASTA other than the chr22 haplotype, and a record of them.

Every record must be uppercase ACGT only, as in the Fig Ops query workloads. With
pair_reverse_complements, records come in pairs, a sequence and then its reverse
complement (the Fig Ops mtb152 queries); each pair is checked and only its first record
is kept, since template strands already cover both strands.
"""
    path = Path(path).resolve(strict=True)
    records = list(fasta(path))
    for name, seq in records:
        if not seq or re.search('[^ACGT]', seq):
            raise ValueError('Source record is empty or not ACGT-only: ' + name)
    if pair_reverse_complements:
        if len(records) % 2:
            raise ValueError('Source records do not come in pairs')
        for (a, x), (b, y) in zip(records[::2], records[1::2]):
            if y != reverse_complement(x):
                raise ValueError(f'{b} is not the reverse complement of {a}')
        segments = records[::2]
    else:
        segments = records
    record = dict(path=str(path), sha256=sha256(path), note=note, records=len(records),
                  pair_reverse_complements=pair_reverse_complements, segments=[n for n, _ in segments],
                  acgt_bases=sum(len(s) for _, s in segments),
                  rule=('ACGT-only records; each forward record of a forward and reverse-complement pair'
                        if pair_reverse_complements else 'every ACGT-only record'))
    return record, segments


def real_lengths(name, real=None):
    """Read lengths of the real workload the set mimics, from its verified offsets."""
    real = Path(real or REAL[name])
    manifest = json.loads((real/'manifest.json').read_text())
    verify(real/'offsets.u64', manifest['files']['offsets.u64']['sha256'], 'real offsets')
    offsets = read_u64(real/'offsets.u64')
    lengths = [offsets[i + 1] - offsets[i] for i in range(len(offsets) - 1)]
    if len(lengths) != manifest['reads'] or offsets[-1] != manifest['bases']:
        raise ValueError('Real workload dimensions differ')
    return dict(path=str(real), manifest_sha256=sha256(real/'manifest.json'), cohort=manifest['cohort'],
                offsets_sha256=manifest['files']['offsets.u64']['sha256']), lengths


def templates(segments, lengths, total, rng):
    """Templates until total bases: length drawn from lengths, start uniform over the positions where it fits, strand 1/2.

Returns [(name, sequence, origin)].
"""
    out, bases = [], 0
    sequences = dict(segments)
    while bases < total:
        length = rng.choice(lengths)
        fits = [(name, len(seq) - length + 1) for name, seq in segments if len(seq) >= length]
        if not fits:
            raise ValueError(f'No ACGT segment holds a template of {length} bases')
        cumulative, running = [], 0
        for _, count in fits:
            running += count
            cumulative.append(running)
        pick = rng.randrange(running)
        k = bisect.bisect_right(cumulative, pick)
        start = pick - (cumulative[k - 1] if k else 0)
        name = fits[k][0]
        sequence = sequences[name][start:start + length]
        strand = '+' if rng.random() < 0.5 else '-'
        if strand == '-':
            sequence = reverse_complement(sequence)
        out.append((f't{len(out)}', sequence, dict(segment=name, begin=start, end=start + length, strand=strand)))
        bases += length
    return out


def parse_art_aln(path):
    """{read id: origin} from an ART ALN file: '>REF\\tREADID\\tPOS\\tSTRAND' headers."""
    out = {}
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                ref, read, position, strand = line[1:].rstrip('\n').split('\t')
                out[read] = dict(segment=ref, begin=int(position), strand=strand)
    return out


def run_tool(argv, work, label, env=None, timeout=TIMEOUT):
    """Run one simulator command in work; stderr and timing are kept as evidence."""
    started = time.time()
    with (work/(label + '.stdout')).open('w') as out, (work/(label + '.stderr')).open('w') as err:
        result = subprocess.run([str(a) for a in argv], cwd=work, stdout=out, stderr=err, timeout=timeout,
                                env=dict(os.environ, **(env or {})))
    row = dict(argv=[str(a) for a in argv], cwd=str(work), returncode=result.returncode,
               wall_seconds=time.time() - started)
    if result.returncode:
        raise RuntimeError(f'{label} failed with {result.returncode}: ' + (work/(label + '.stderr')).read_text()[-2000:])
    return row


def simulate_short(segments_path, segments, work, seed, oversample, acgt_bases, art=None):
    art = art or binary('art_illumina')
    fold = round(TARGET * oversample / acgt_bases, 4)
    argv = [art['path'], '-ss', 'HS25', '-i', segments_path, '-l', SHORT_LENGTH, '-f', fold, '-rs', seed,
            '-o', work/'art']
    command = run_tool(argv, work, 'art')
    origins = parse_art_aln(work/'art.aln')
    reads = []
    for name, seq, _ in fastq(work/'art.fq'):
        origin = origins[name]
        reads.append((seq, dict(origin, end=origin['begin'] + SHORT_LENGTH, source=name)))
    return reads, dict(tools=[art], commands=[command], fold_coverage=fold,
                       profile='ART built-in HiSeq 2500 profile (-ss HS25; ART reports "HiSeq 2500 Length 126 R1"), '
                               'single-end, default indel rates, strand of each read chosen by ART with probability 1/2')


def ccs_accuracy(quals):
    """Predicted accuracy of ccs reads from their base qualities (Phred+33), per read and over all bases.

A read's predicted accuracy is one minus its mean per-base error probability, the definition
of the ccs rq value that the FASTQ output does not carry.
"""
    error = [10 ** (-(q - 33) / 10) for q in range(128)]
    per_read, errors, bases = [], 0.0, 0
    for qual in quals:
        e = sum(error[ord(c)] * n for c, n in collections.Counter(qual).items())
        per_read.append(1 - e / len(qual))
        errors += e
        bases += len(qual)
    if not per_read:
        return {}
    ranked = sorted(per_read)
    pick = lambda f: ranked[min(len(ranked) - 1, int(f * len(ranked)))]
    return dict(reads=len(per_read), mean_read=sum(per_read) / len(per_read), base_weighted=1 - errors / bases,
                q01=pick(.01), q10=pick(.1), median=pick(.5), min=ranked[0],
                definition='1 - mean per-base error probability from ccs base qualities (Phred+33), as ccs rq')


def simulate_pacbio(segments, work, seed, oversample, lengths, threads, pbsim=None, ccs=None, samtools=None,
                    accuracy=False):
    """With accuracy, details also record pbsim3's subread statistics and the ccs reads' predicted accuracy."""
    pbsim, ccs, samtools = pbsim or binary('pbsim'), ccs or binary('ccs'), samtools or binary('samtools')
    rng = random.Random(seed)
    made = templates(segments, lengths, TARGET * oversample, rng)
    with (work/'templates.fa').open('w') as f:
        for name, seq, _ in made:
            f.write(f'>{name}\n{seq}\n')
    # ccs names its per-thread temporary BAM files thread.N_0.ccs.fastq.bam in TMPDIR, so two jobs on
    # one node that share a TMPDIR overwrite each other's files; the job's own work directory avoids that.
    env = dict(PATH=str(Path(samtools['path']).parent) + os.pathsep + os.environ.get('PATH', ''), TMPDIR=str(work))
    model = MODEL['qshmm-rsii']
    one = run_tool([pbsim['path'], '--strategy', 'templ', '--template', work/'templates.fa', '--method', 'qshmm',
                    '--qshmm', model, '--pass-num', HIFI_PASSES, '--seed', seed, '--prefix', work/'pbsim'],
                   work, 'pbsim', env)
    two = run_tool([ccs['path'], '-j', threads, work/'pbsim.bam', work/'ccs.fastq.gz'], work, 'ccs', env)
    reads, quals = [], []
    for name, seq, qual in fastq(work/'ccs.fastq.gz'):
        if accuracy:
            quals.append(qual)
        parts = name.split('/')
        if len(parts) != 3 or parts[2] != 'ccs':
            raise ValueError('Unexpected ccs read name ' + name)
        template = made[int(parts[1]) - 1]
        reads.append((seq, dict(template[2], source=name, template=template[0], template_strand=template[2]['strand'])))
    extra = dict(accuracy=dict(subreads=pbsim_stats(work/'pbsim.stderr'), ccs_predicted=ccs_accuracy(quals))) if accuracy else {}
    return reads, dict(tools=[pbsim, ccs, samtools], commands=[one, two], model=dict(path=str(model), sha256=sha256(model)),
                       **extra,
                       templates=dict(count=len(made), bases=sum(len(s) for _, s, _ in made), seed=seed,
                                      sha256=sha256(work/'templates.fa'),
                                      rule='length drawn with replacement from the real PacBio workload, start uniform '
                                           'over the ACGT-segment positions where it fits, strand + or - with probability 1/2'),
                       ccs_report=(work/'ccs.ccs_report.txt').read_text() if (work/'ccs.ccs_report.txt').exists() else None,
                       profile=f'pbsim3 template strategy, QSHMM-RSII.model with its default pass accuracy (0.85) and '
                               f'difference ratio, --pass-num {HIFI_PASSES}; ccs 6.4.0 defaults (min passes 3, min rq 0.99). '
                               f'The read strand is the template strand, or its reverse if ccs reports the other orientation.')


def pbsim_stats(path):
    """Read count, length and accuracy lines of a pbsim3 stderr report, from its last statistics block."""
    out = {}
    for line in Path(path).read_text().splitlines():
        key, sep, value = line.partition(':')
        key = key.strip().rstrip('.').strip()
        if sep and key in ('read num', 'read length mean (SD)', 'read length min', 'read length max',
                           'read accuracy mean (SD)', 'substitution rate', 'insertion rate', 'deletion rate'):
            out[key] = value.strip()
    return out


def parse_templ_maf(path):
    """{read name: template name} from a pbsim3 template-strategy MAF."""
    out, ref = {}, None
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith('a'):
                ref = None
            elif line.startswith('s '):
                name = line.split()[1]
                if ref is None:
                    ref = name
                else:
                    out[name] = ref
    return out


def simulate_ont(segments, work, seed, oversample, lengths, pbsim=None):
    pbsim = pbsim or binary('pbsim')
    made = templates(segments, lengths, TARGET * oversample, random.Random(seed))
    with (work/'templates.fa').open('w') as f:
        for name, seq, _ in made:
            f.write(f'>{name}\n{seq}\n')
    model = MODEL['errhmm-ont-hq']
    command = run_tool([pbsim['path'], '--strategy', 'templ', '--template', work/'templates.fa', '--method', 'errhmm',
                        '--errhmm', model, '--accuracy-mean', ONT_ACCURACY_MEAN, '--seed', seed, '--prefix', work/'pbsim'], work, 'pbsim')
    origin = {name: o for name, _, o in made}
    template_of = parse_templ_maf(work/'pbsim.maf.gz')
    reads = []
    for read, seq, _ in fastq(work/'pbsim.fq.gz'):
        template = template_of[read]
        reads.append((seq, dict(origin[template], source=read, template=template)))
    if len(reads) != len(made):
        raise ValueError('pbsim3 did not return one read per template')
    return reads, dict(tools=[pbsim], commands=[command], model=dict(path=str(model), sha256=sha256(model)),
                       accuracy=pbsim_stats(work/'pbsim.stderr'),
                       templates=dict(count=len(made), bases=sum(len(s) for _, s, _ in made), seed=seed,
                                      sha256=sha256(work/'templates.fa'),
                                      rule='length drawn with replacement from the real ONT workload, start uniform '
                                           'over the ACGT-segment positions where it fits, strand + or - with probability 1/2'),
                       profile='pbsim3 template strategy, --method errhmm with ERRHMM-ONT-HQ.model (its standard ONT '
                               'model, recommended for read sets of 90 percent accuracy or more), --accuracy-mean '
                               f'{ONT_ACCURACY_MEAN}, the model\'s ceiling; one read per template, on the template\'s strand')


def finish(name, reads, destination, seed, details, haplotype_manifest, segments_table=None, source=None, cohort=None):
    """Drop reads with non-ACGT, shuffle, take whole reads to TARGET, renumber and write the read set.

With source (a load_source record), the manifest names that source in place of the chr22
haplotype; cohort replaces the set's chr22 cohort.
"""
    destination = Path(destination)
    kept, dropped = [], dict(reads=0, bases=0)
    for seq, origin in reads:
        upper = seq.upper()
        if not set(upper) <= ACGT:
            dropped['reads'] += 1
            dropped['bases'] += len(seq)
            continue
        kept.append((upper, origin))
    random.Random(seed).shuffle(kept)
    chosen, bases = [], 0
    for row in kept:
        if bases >= TARGET:
            break
        chosen.append(row)
        bases += len(row[0])
    if bases < TARGET:
        raise ValueError(f'Only {bases} simulated ACGT bases; the target is {TARGET}')
    destination.mkdir(parents=True, exist_ok=False)
    offsets = array.array('Q', [0])
    strands = {}
    with (destination/'reads.fa').open('x') as f, (destination/'origins.tsv').open('x') as o:
        o.write('name\tlength\tsegment\tbegin\tend\tstrand\tsource\n')
        for i, (seq, origin) in enumerate(chosen):
            f.write(f'>q{i}\n{seq}\n')
            offsets.append(offsets[-1] + len(seq))
            strands[origin.get('strand')] = strands.get(origin.get('strand'), 0) + 1
            o.write(f"q{i}\t{len(seq)}\t{origin.get('segment', '')}\t{origin.get('begin', '')}\t{origin.get('end', '')}\t"
                    f"{origin.get('strand', '')}\t{origin.get('source', '')}\n")
    with (destination/'one.fa').open('x') as f:
        f.write(f'>q0\n{chosen[0][0]}\n')
    if sys.byteorder != 'little':
        offsets.byteswap()
    (destination/'offsets.u64').write_bytes(offsets.tobytes())
    files = ('reads.fa', 'one.fa', 'offsets.u64', 'origins.tsv')
    if source is None:
        origin_entry = dict(source_haplotype=dict(path=str(HAPLOTYPE), sample=haplotype_manifest['sample'],
                                                  manifest_sha256=haplotype_manifest.get('manifest_sha256'),
                                                  acgt_bases=haplotype_manifest['acgt_bases']))
    else:
        origin_entry = dict(source_segments=source)
    manifest = dict(
        format=FORMAT, status='reads-prepared', set=name, cohort=cohort or SETS[name]['cohort'], platform=SETS[name]['platform'],
        reads=len(chosen), bases=bases, one_bases=len(chosen[0][0]), read_filter=FILTER, target_bases=TARGET,
        lengths=quantiles([len(s) for s, _ in chosen]), strands=strands,
        selection='simulated reads uppercased; reads holding any non-ACGT base dropped; the rest shuffled with the '
                  'seed and whole reads taken in that order until at least the target; renumbered q0, q1, ...',
        simulated=dict(reads=len(reads), kept_acgt=len(kept), dropped_non_acgt=dropped),
        seed=seed, **origin_entry, simulation=details, host=socket.gethostname(), job_id=os.environ.get('SLURM_JOB_ID'), created=time.time(),
        files={n: dict(bytes=(destination/n).stat().st_size, sha256=sha256(destination/n)) for n in files})
    save_new(destination/'manifest.json', manifest)
    return manifest


def simulate(name, threads=8, destination=None, work=None, source=None, seed=None, cohort=None, real_workload=None):
    """Simulate one read set from the chr22 haplotype or, with source, from (record, segments) of load_source.

A source run needs a destination and may set its own seed and cohort; the recipe (simulator,
model, template lengths, oversampling, filter and target) is the set's.
"""
    if source is not None and destination is None:
        raise ValueError('A simulation from another source needs a destination')
    destination = Path(destination or reads_dir(name))
    if destination.exists():
        raise FileExistsError(destination)
    spec = SETS[name]
    seed = spec['seed'] if seed is None else seed
    if source is None:
        hap, segments = load_haplotype()
        hap['manifest_sha256'] = sha256(HAPLOTYPE/'manifest.json')
        label = name
    else:
        hap, segments = source
        label = destination.name
    real, lengths = real_lengths(name, real_workload)
    if name == 'short' and set(lengths) != {SHORT_LENGTH}:
        raise ValueError('The short-read workload is not all 125-base reads')
    work = Path(work or Path(os.environ.get('TMPDIR', '/tmp'))/f"msbakery-1012-chr22sim-{label}-{os.environ.get('SLURM_JOB_ID', os.getpid())}")
    work.mkdir(parents=True)
    evidence = P/'chr22-sim-jobs'/f"{label}-{os.environ.get('SLURM_JOB_ID', 'local')}"
    evidence.mkdir(parents=True)
    try:
        segments_path = HAPLOTYPE/'acgt.fa'
        if name == 'short' and source is not None:
            segments_path = work/'segments.fa'
            with segments_path.open('x') as f:
                for n, s in segments:
                    f.write(f'>{n}\n{s}\n')
        if name == 'short':
            reads, details = simulate_short(segments_path, segments, work, seed, spec['oversample'], hap['acgt_bases'])
        elif name == 'pacbio':
            # The chr22 manifests predate the accuracy record, so only other sources add it.
            reads, details = simulate_pacbio(segments, work, seed, spec['oversample'], lengths, threads,
                                             accuracy=source is not None)
        else:
            reads, details = simulate_ont(segments, work, seed, spec['oversample'], lengths)
        details.update(real_workload=real, real_lengths=quantiles(lengths), oversample=spec['oversample'])
        if source is None:
            manifest = finish(name, reads, destination, seed, details, hap)
        else:
            manifest = finish(name, reads, destination, seed, details, None, source=hap, cohort=cohort)
        return manifest
    finally:
        for path in work.glob('*'):
            if path.suffix in ('.stderr', '.stdout', '.txt'):
                shutil.copyfile(path, evidence/path.name)
        shutil.rmtree(work, ignore_errors=True)


# Historical chr22 cells: their staging.json files hashed the binaries, FMDs and legacy indexes.
CELLS = ROOT/'results/real-hprc2-1d580b6/chr22/cells'


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(16 << 20), b''):
            h.update(block)
    return h.hexdigest()


def staged_record(batch, member, suffix):
    """(path, bytes, sha256) of the one file a batch member staged whose source ends with suffix."""
    rows = [f for f in json.loads((Path(batch)/member/'staging.json').read_text())['files']
            if str(f.get('source', '')).endswith(suffix)]
    if len(rows) != 1 or len(rows[0].get('sha256', '')) != 64:
        raise ValueError(f'No single verified staging record for {suffix} in {member}')
    return dict(path=rows[0]['source'], bytes=rows[0]['bytes'], sha256=rows[0]['sha256'])


def rope_records(reference):
    """The ropebwt3 binary and FMD the chr22 runs stage, from the historical staging records."""
    batch = CELLS/reference/'human_short'
    rope = staged_record(batch, 'ropebwt3', '/bin/ropebwt3')
    fmd = staged_record(batch, 'ropebwt3', '/ropebwt3/index.fmd')
    if rope['sha256'] != ROPE_SHA256:
        raise ValueError('ropebwt3 differs from the frozen comparator')
    return rope, fmd


def bounded(truth, offsets):
    """Every ACGT base occurs in the reference and no match runs past its read."""
    with Path(truth).open('rb') as f:
        for i in range(len(offsets) - 1):
            n = offsets[i + 1] - offsets[i]
            ms = array.array('I')
            ms.frombytes(f.read(4 * n))
            if sys.byteorder != 'little':
                ms.byteswap()
            if len(ms) != n or any(not 1 <= v <= n - j for j, v in enumerate(ms)):
                raise ValueError(f'MS truth out of range for q{i}')


def truth(name, reference, source=None, destination=None, work=None, records=rope_records, timeout=6 * 3600):
    source = Path(source or reads_dir(name))
    destination = Path(destination or workload_dir(reference, name))
    if destination.exists():
        raise FileExistsError(destination)
    reads = json.loads((source/'manifest.json').read_text())
    if reads.get('status') != 'reads-prepared' or reads.get('set') != name or reads.get('read_filter') != FILTER:
        raise ValueError('Prepared simulated reads required')
    for n, record in reads['files'].items():
        verify(source/n, record['sha256'], n)
    rope, fmd = records(reference)
    verify(rope['path'], rope['sha256'], 'ropebwt3')
    verify(fmd['path'], fmd['sha256'], 'FMD')
    offsets = read_u64(source/'offsets.u64')
    if len(offsets) != reads['reads'] + 1 or offsets[-1] != reads['bases']:
        raise ValueError('Read boundaries differ from the manifest')
    work = Path(work or P/'truth-jobs'/f"chr22sim-{reference}-{name}-{os.environ.get('SLURM_JOB_ID', 'local')}")
    work.mkdir(parents=True, exist_ok=False)
    started = time.time()
    # The historical chr22 truth command (job 31647376): one thread, every SMEM of length 1 or more.
    row = common.execute([rope['path'], 'mem', '-t1', '-l1', fmd['path'], source/'reads.fa'], work, 'truth', timeout)
    partial = work/'truth.u32'
    row['check'] = check('ropebwt3', work, row, partial, offsets, write=True)
    bounded(partial, offsets)
    Path(row['stdout']).unlink()
    row['raw_output_retained'] = False
    common.save(work/'truth.command.json', row)
    destination.mkdir(parents=True)
    for n in ('reads.fa', 'one.fa', 'offsets.u64'):
        shutil.copyfile(source/n, destination/n)
        verify(destination/n, reads['files'][n]['sha256'], 'copied ' + n)
    with partial.open('rb') as src, (destination/'truth.u32').open('xb') as dst:
        shutil.copyfileobj(src, dst, 16 << 20)
        dst.flush()
        os.fsync(dst.fileno())
    digest = sha256(destination/'truth.u32')
    if digest != sha256(partial):
        raise ValueError('Truth copy differs')
    partial.unlink()
    verdict = dict(passed=True, command=row, reference=reference, reads_manifest_sha256=sha256(source/'manifest.json'),
                   sha256=digest, ropebwt3=rope, fmd=fmd,
                   derivation='ropebwt3 mem -t1 -l1 SMEMs converted to per-position MS by ms2exp.validate.ms.check; '
                              'every value checked to satisfy 1 <= ms[j] <= length - j',
                   job_id=os.environ.get('SLURM_JOB_ID'), started=started, finished=time.time())
    save_new(destination/'truth.json', verdict)
    names = ('reads.fa', 'one.fa', 'offsets.u64', 'truth.u32', 'truth.json')
    manifest = dict(status='prepared', reference=reference, cohort=reads['cohort'], set=name, reads=reads['reads'],
                    bases=reads['bases'], one_bases=reads['one_bases'], read_filter=FILTER,
                    selection=f"simulated {reads['platform']} reads of {source} (same records, names and order)",
                    reads_workload=str(source), reads_manifest_sha256=sha256(source/'manifest.json'),
                    truth=(f"ropebwt3 -t1 -l1 (binary {rope['sha256'][:12]}) on the {reference} FMD "
                           f"({fmd['sha256'][:12]}), job {os.environ.get('SLURM_JOB_ID')}"),
                    truth_metadata_sha256=sha256(destination/'truth.json'), reuse_old_timings=False,
                    files={n: dict(bytes=(destination/n).stat().st_size, sha256=sha256(destination/n)) for n in names})
    save_new(destination/'manifest.json', manifest)
    return manifest
