#!/bin/bash
# Prepare the toy end-to-end run of the freeze workflow: tools built from the ms2 commit, the toy's
# ms2 indexes, its two cut read sets, the inputs catalog, the toy manifest and its inputs tree.
# Usage: toy_prep.sh TOY_DIR [--data DIR] [--ms2-repo URL] [--ms2-commit SHA] [ms2exp.build.tools options ...]
# Every input comes from the data directory (--data, else MS2EXP_DATA), which
# `python -m ms2exp.fetch --set toy` fills; this script checks it first. The ms2 repository and
# commit default to tools.lock.json (langmead-lab/ms2 at paper-v1); the repository may also be a
# local path or a git bundle (toy/make_source_bundles.sh), and the options after it go to
# ms2exp.build.tools (for example --submodule-repo NAME URL). THREADS (default: the core count, at
# most 16) sets the threads of the index and SSA builds.
# The toy builds the mtb152 indexes from their precursors with the ms2 build commands of the
# all-AWS build jobs (ms2exp.manifest build_matrix): min with --minima, psi, and phi and full with
# --phi and --minima as for an FMD, plus ftabs for min and psi and a full index with seqids, and
# the chr22-h10 min index the legacy job's ms2 cell reads. It also makes what locate needs: a
# ropebwt3 .len.gz from the FMD itself and an s8 SSA. It cuts two read sets with their truth
# sliced to match: 2 Mbases of the mtb152 PacBio set for locate, and 4 Mbases of the chr22-h10
# simulated short reads for the legacy job, which runs MONI, PHONI and AugPHONI on their chr22-h10
# indexes with the pinned binaries (D1).
set -euo pipefail
T=$(mkdir -p "$1" && cd "$1" && pwd -P)
shift
REPO=$(dirname "$(dirname "$(readlink -f "$0")")")
PY="env PYTHONPATH=$REPO${PYTHONPATH:+:$PYTHONPATH} python3"
lock() { python3 -c "import json, sys; print(json.load(open('$REPO/tools.lock.json'))['ms2'][sys.argv[1]])" "$1"; }
DATA=${MS2EXP_DATA:-}
MS2_REPO=$(lock repo)
COMMIT=$(lock commit)
while [ $# -gt 0 ]; do
    case $1 in
        --data) DATA=$2; shift 2 ;;
        --ms2-repo) MS2_REPO=$2; shift 2 ;;
        --ms2-commit) COMMIT=$2; shift 2 ;;
        *) break ;;
    esac
done
[ -n "$DATA" ] || { echo "Give --data DIR or set MS2EXP_DATA (python -m ms2exp.fetch --set toy fills it)" >&2; exit 1; }
DATA=$(cd "$DATA" && pwd -P)
THREADS=${THREADS:-$(n=$(nproc); echo $((n < 16 ? n : 16)))}

$PY -m ms2exp.fetch --data "$DATA" --set toy --verify-only > "$T/fetch-verify.log"
ROPE=$DATA/tools/ropebwt3-cfad8475/ropebwt3
for tool in "$ROPE" "$DATA"/tools/{moni/moni-ms,phoni/phoni,augphoni/aug_phoni}; do
    [ -x "$tool" ] || { echo "$tool is not executable: rerun python3 -m ms2exp.fetch --set toy, which sets its mode" >&2; exit 1; }
done
# The pinned ropebwt3 loads the OpenMP runtime shipped beside it.
RL=$DATA/tools/runtime${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
M=$DATA/precursors/mtb152
C=$DATA/precursors/chr22-h10
FMD=$DATA/ropebwt3/mtb152/index.fmd
R=$DATA/reads
X=$T/extra
mkdir -p "$X"

$PY -m ms2exp.build.tools --lock "$REPO/tools.lock.json" --ms2-repo "$MS2_REPO" --ms2-commit "$COMMIT" \
    --root "$T/built" --ropebwt3 "$ROPE" --runtime "$DATA/tools/runtime" \
    --legacy moni="$DATA/tools/moni/moni-ms" --legacy phoni="$DATA/tools/phoni/phoni" \
    --legacy augphoni="$DATA/tools/augphoni/aug_phoni" --work "$T/ms2-src-${COMMIT:0:7}" "$@"
MS=$T/built/tools/ms2-${COMMIT:0:7}/ms2-bench

# Sequence lengths in FMD order: sequence 2k is the forward strand of contig k.
if [ ! -s "$X/mtb152.fmd.len.gz" ]; then
    n=$(env LD_LIBRARY_PATH="$RL" "$ROPE" stat "$FMD" | awk '$1 == "sequences" || /Number of sequences/ {print $NF}' | head -1)
    n=${n:-288}
    seq 0 2 $((n - 2)) | xargs -P "$THREADS" -I{} env LD_LIBRARY_PATH="$RL" sh -c \
        "'$ROPE' get '$FMD' {} | awk -v i={} 'NR > 1 {l += length(\$0)} END {print i/2 \"\t\" l}'" \
        | sort -n | awk '{print "s" $1 "\t" $2}' | gzip -n > "$X/mtb152.fmd.len.gz"
fi
[ -s "$X/mtb152.fmd.s8.ssa" ] || env LD_LIBRARY_PATH="$RL" "$ROPE" ssa -o "$X/mtb152.fmd.s8.ssa" -s8 -t"$THREADS" "$FMD"
# ms2 indexes at this commit, built as build_matrix's index_steps builds them.
I=$X/ms2-${COMMIT:0:7}
mkdir -p "$I"
[ -s "$I/min.idx" ] || "$MS" build "$M/bwt.heads" "$M/bwt.len" "$I/min.idx" --layout min --minima "$M/min"
[ -s "$I/psi.idx" ] || "$MS" build "$M/bwt.heads" "$M/bwt.len" "$I/psi.idx" --layout psi
for layout in phi full; do
    [ -s "$I/$layout.idx" ] || "$MS" build "$M/bwt.heads" "$M/bwt.len" "$I/$layout.idx" --layout $layout \
        --phi "$M/phi" --minima "$M/min"
done
for layout in min psi; do
    [ -s "$I/$layout.idx.ftab" ] || "$MS" ftab-build "$I/$layout.idx" "$I/$layout.idx.ftab" --auto -p "$THREADS"
done
[ -s "$I/full.seqids.idx" ] || "$MS" index-add-seq-ids "$I/full.idx" "$X/mtb152.fmd.len.gz" "$I/full.seqids.idx"
# The chr22-h10 min index (about 20 s and 0.6 GB).
[ -s "$I/chr22-h10.min.idx" ] || "$MS" build "$C/bwt.heads" "$C/bwt.len" "$I/chr22-h10.min.idx" --layout min --minima "$C/min"

# cut_set SRC DST BASES: the first reads of SRC up to BASES bases, with their slice of the MS
# truth; one.fa (the first read, for the legacy one-read run) is copied when SRC has it.
cut_set() {
    mkdir -p "$2"
    python3 - "$1" "$2" "$3" <<'EOF'
import array, shutil, sys
from pathlib import Path
src, dst, limit = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
bounds = array.array('Q'); bounds.frombytes((src/'offsets.u64').read_bytes())
n = next(i for i in range(len(bounds)) if bounds[i] >= limit)
with (src/'reads.fa').open() as f, (dst/'reads.fa').open('w') as out:
    seen = 0
    for line in f:
        if line.startswith('>'):
            if seen == n:
                break
            seen += 1
        out.write(line)
array.array('Q', bounds[:n + 1]).tofile((dst/'offsets.u64').open('wb'))
with (src/'truth.u32').open('rb') as f:
    (dst/'truth.u32').write_bytes(f.read(bounds[n] * 4))
if (src/'one.fa').exists():
    shutil.copyfile(src/'one.fa', dst/'one.fa')
    first = (dst/'reads.fa').read_text().split('>', 2)[1]
    if (dst/'one.fa').read_text() != '>' + first:
        raise SystemExit('one.fa is not the first read of ' + str(src))
print(f'{dst.name}: {n} reads, {bounds[n]} bases')
EOF
}
# Locate reads: the first PacBio reads up to 2 Mbases.
L=$X/mtb152-locate-2mb-v2
cut_set "$R/mtb152/mtb152-pacbio-60mb-v2" "$L" 2000000
# Legacy reads: the first chr22-h10 simulated short reads up to 4 Mbases.
S=$X/chr22-h10-simshort-4mb-v2
cut_set "$R/chr22-h10/simshort-60mb-v2" "$S" 4000000

# The fetched inputs the toy reads directly (not the tools, which ms2exp.build.tools installed
# under built/tools, nor the precursors, which only the builds above read), then what it made.
python3 - "$DATA/locations.json" "$T/toy-locations.json" <<EOF
import json, sys
fetched = json.load(open(sys.argv[1]))['locations']
locations = {n: p for n, p in fetched.items() if not n.startswith(('tools/', 'precursors/'))}
for f in ('reads.fa', 'truth.u32', 'offsets.u64'):
    locations[f'reads/mtb152/mtb152-locate-2mb-v2/{f}'] = f'$L/{f}'
for f in ('reads.fa', 'truth.u32', 'offsets.u64', 'one.fa'):
    locations[f'reads/chr22-h10/simshort-4mb-v2/{f}'] = f'$S/{f}'
locations.update({
    'tools': '$T/built/tools',
    'ms2/mtb152/index.ms.orbitidx': '$I/min.idx',
    'ms2/mtb152/index.ms.orbitidx.ftab': '$I/min.idx.ftab',
    'ms2/mtb152/index.tms.psi.orbitidx': '$I/psi.idx',
    'ms2/mtb152/index.tms.psi.orbitidx.ftab': '$I/psi.idx.ftab',
    'ms2/mtb152/index.tms.phi.orbitidx': '$I/phi.idx',
    'ms2/mtb152/index.tms.full.orbitidx': '$I/full.idx',
    'ms2/mtb152/index.tms.full.seqids.orbitidx': '$I/full.seqids.idx',
    'ms2/chr22-h10/index.ms.orbitidx': '$I/chr22-h10.min.idx',
    'ropebwt3/mtb152/index.fmd.len.gz': '$X/mtb152.fmd.len.gz',
    'ropebwt3/mtb152/index.fmd.s8.ssa': '$X/mtb152.fmd.s8.ssa',
})
json.dump(dict(locations=dict(sorted(locations.items()))), open(sys.argv[2], 'w'), indent=1)
EOF
$PY -m ms2exp.catalog hash --locations "$T/toy-locations.json" --out "$T/toy-inputs.json"
$PY -m ms2exp.catalog link --locations "$T/toy-locations.json" --root "$T/inputs-root"
$PY -m ms2exp.manifest --matrix toy --inputs "$T/toy-inputs.json" --ms2-repo "$MS2_REPO" --ms2-commit "$COMMIT" \
    --note "toy end-to-end run; mtb152 and chr22-h10 indexes built at the ms2 commit" -o "$T/manifest-toy.json"
