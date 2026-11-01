#!/bin/bash
# Prepare the toy end-to-end run of the freeze workflow on Rockfish (devlangmead1): mtb152 inputs,
# tools built from the ms2 commit, the inputs catalog, the toy manifest and its inputs tree.
# Usage: toy_prep.sh TOY_DIR MS2_REPO MS2_COMMIT [--submodule-repo NAME URL ...]
# MS2_REPO is a URL, a local repository path or a git bundle; the options after MS2_COMMIT go to
# ms2exp.build.tools. Inputs are Rockfish paths (lab-internal).
# The toy builds the mtb152 indexes from their precursors with the ms2 build commands of the
# all-AWS build jobs (ms2exp.manifest build_matrix): min with --minima, psi, and phi and full with
# --phi and --minima as for an FMD, plus ftabs for min and psi and a full index with seqids. It
# also makes what locate needs: a ropebwt3 .len.gz from the FMD itself, an s8 SSA, and a 2 Mb
# locate read set cut from the mtb152 PacBio set with its truth sliced to match. The legacy job
# reuses the chr22-h10 MONI, PHONI and AugPHONI indexes and their pinned binaries (D1), the
# 1d580b6 chr22-h10 min index (which ms2 reads), and a 4 Mb cut of the chr22-h10 simulated
# short-read set with its truth and one-read file.
set -euo pipefail
T=$(readlink -f "$1")
MS2_REPO=$2
COMMIT=$3
shift 3
REPO=$(dirname "$(dirname "$(readlink -f "$0")")")
PY="env PYTHONPATH=$REPO${PYTHONPATH:+:$PYTHONPATH} python3"
P=/scratch16/blangme2/ms-bakery/results/campaign-20261012
ROPE=/scratch16/blangme2/ms-bakery/results/real-hprc2-1d580b6/bin/ropebwt3
D=/vast/blangme2/terams-fig/build/mtb152/orbit
FMD=/vast/blangme2/nbrown/ms_lcp/data/mtb152/mtb152.95%.fmd
X=$T/extra
mkdir -p $X
# chr22-h10 legacy indexes and binaries, as staged by the real-hprc2 chr22 cells.
LADDER=/scratch16/blangme2/ms-bakery/results/ladder/chr22-h10
LEGBIN=/scratch16/blangme2/ms-bakery/results/real-hprc2-1d580b6/chr22/bin
CHR22_MS=/scratch16/blangme2/ms-bakery/results/orbit-1d580b6-chr22/chr22-h10/Orbit-ms/index.orbitidx
SHORT=$P/inputs/chr22_sim/chr22-h10-short-60mb-v1
# Pins of chr22-legacy-binary-verification.json.
sha256sum -c --quiet <<SUMS
4ea254ddf5921341a45029fe94279618db38f7db47c26b57db9695ea586efb3f  $LEGBIN/moni-ms
aeb7fb3b4ec4f82203f015102fe3429345f2ec166ee7ce9d5a2644bc23e83c58  $LEGBIN/phoni
79d50bf5bf4997ec35587eb1b0c17bdce7c728ddfd472305272c73a7317e862b  $LEGBIN/aug_phoni
SUMS

$PY -m ms2exp.build.tools --ms2-repo $MS2_REPO --ms2-commit $COMMIT --root $T/built --ropebwt3 $ROPE \
    --legacy moni=$LEGBIN/moni-ms --legacy phoni=$LEGBIN/phoni --legacy augphoni=$LEGBIN/aug_phoni \
    --work $T/ms2-src-${COMMIT:0:7} "$@"
MS=$T/built/tools/ms2-${COMMIT:0:7}/ms2-bench

# Sequence lengths in FMD order: sequence 2k is the forward strand of contig k.
if [ ! -s $X/mtb152.fmd.len.gz ]; then
    n=$($ROPE stat $FMD 2>/dev/null | awk '$1 == "sequences" || /Number of sequences/ {print $NF}' | head -1)
    n=${n:-288}
    seq 0 2 $((n - 2)) | xargs -P 16 -I{} sh -c "$ROPE get $FMD {} | awk -v i={} 'NR > 1 {l += length(\$0)} END {print i/2 \"\t\" l}'" \
        | sort -n | awk '{print "s" $1 "\t" $2}' | gzip -n > $X/mtb152.fmd.len.gz
fi
[ -s $X/mtb152.fmd.s8.ssa ] || $ROPE ssa -o $X/mtb152.fmd.s8.ssa -s8 -t16 $FMD
# mtb152 indexes at this commit, built as build_matrix's index_steps builds them.
I=$X/ms2-${COMMIT:0:7}
mkdir -p $I
[ -s $I/min.idx ] || $MS build $D/mtb152.bwt.heads $D/mtb152.bwt.len $I/min.idx --layout min --minima $D/mtb152.min
[ -s $I/psi.idx ] || $MS build $D/mtb152.bwt.heads $D/mtb152.bwt.len $I/psi.idx --layout psi
for layout in phi full; do
    [ -s $I/$layout.idx ] || $MS build $D/mtb152.bwt.heads $D/mtb152.bwt.len $I/$layout.idx --layout $layout \
        --phi $D/mtb152.phi --minima $D/mtb152.min
done
for layout in min psi; do
    [ -s $I/$layout.idx.ftab ] || $MS ftab-build $I/$layout.idx $I/$layout.idx.ftab --auto -p 16
done
[ -s $I/full.seqids.idx ] || $MS index-add-seq-ids $I/full.idx $X/mtb152.fmd.len.gz $I/full.seqids.idx

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
L=$X/mtb152-locate-2mb-v1
cut_set $P/inputs/fig_ops2/mtb152-pacbio-60mb-v1 $L 2000000
# Legacy reads: the first chr22-h10 simulated short reads up to 4 Mbases.
S=$X/chr22-h10-simshort-4mb-v1
cut_set $SHORT $S 4000000

cat > $T/toy-locations.json <<EOF
{"locations": {
  "tools": "$T/built/tools",
  "reads/mtb152/mtb152-heldout16-p100k-v1/reads.fa": "$P/inputs/fig_ops/mtb152-heldout16-p100k-v1/reads.fa",
  "reads/mtb152/mtb152-heldout16-p100k-v1/truth.u32": "$P/inputs/fig_ops/mtb152-heldout16-p100k-v1/truth.u32",
  "reads/mtb152/mtb152-heldout16-p100k-v1/offsets.u64": "$P/inputs/fig_ops/mtb152-heldout16-p100k-v1/offsets.u64",
  "reads/mtb152/mtb152-pacbio-60mb-v1/reads.fa": "$P/inputs/fig_ops2/mtb152-pacbio-60mb-v1/reads.fa",
  "reads/mtb152/mtb152-pacbio-60mb-v1/truth.u32": "$P/inputs/fig_ops2/mtb152-pacbio-60mb-v1/truth.u32",
  "reads/mtb152/mtb152-pacbio-60mb-v1/offsets.u64": "$P/inputs/fig_ops2/mtb152-pacbio-60mb-v1/offsets.u64",
  "reads/mtb152/mtb152-locate-2mb-v1/reads.fa": "$L/reads.fa",
  "reads/mtb152/mtb152-locate-2mb-v1/truth.u32": "$L/truth.u32",
  "reads/mtb152/mtb152-locate-2mb-v1/offsets.u64": "$L/offsets.u64",
  "ms2/mtb152/index.ms.orbitidx": "$I/min.idx",
  "ms2/mtb152/index.ms.orbitidx.ftab": "$I/min.idx.ftab",
  "ms2/mtb152/index.tms.psi.orbitidx": "$I/psi.idx",
  "ms2/mtb152/index.tms.psi.orbitidx.ftab": "$I/psi.idx.ftab",
  "ms2/mtb152/index.tms.phi.orbitidx": "$I/phi.idx",
  "ms2/mtb152/index.tms.full.orbitidx": "$I/full.idx",
  "ms2/mtb152/index.tms.full.seqids.orbitidx": "$I/full.seqids.idx",
  "ropebwt3/mtb152/index.fmd": "$FMD",
  "ropebwt3/mtb152/index.fmd.len.gz": "$X/mtb152.fmd.len.gz",
  "ropebwt3/mtb152/index.fmd.s8.ssa": "$X/mtb152.fmd.s8.ssa",
  "reads/chr22-h10/simshort-4mb-v1/reads.fa": "$S/reads.fa",
  "reads/chr22-h10/simshort-4mb-v1/truth.u32": "$S/truth.u32",
  "reads/chr22-h10/simshort-4mb-v1/offsets.u64": "$S/offsets.u64",
  "reads/chr22-h10/simshort-4mb-v1/one.fa": "$S/one.fa",
  "ms2/chr22-h10/index.ms.orbitidx": "$CHR22_MS",
  "moni/chr22-h10/index.thrbv.ms": "$LADDER/moni/index.thrbv.ms",
  "moni/chr22-h10/index.plain.slp": "$LADDER/moni/index.plain.slp",
  "moni/chr22-h10/index.idx": "$LADDER/moni/index.idx",
  "phoni/chr22-h10/reference.txt.phoni": "$LADDER/phoni/reference.txt.phoni",
  "phoni/chr22-h10/reference.txt.plain.slp": "$LADDER/phoni/reference.txt.plain.slp",
  "phoni/chr22-h10/reference.txt": "$LADDER/phoni/reference.txt",
  "augphoni/chr22-h10/reference.txt.aug": "$LADDER/aug_phoni/reference.txt.aug",
  "augphoni/chr22-h10/reference.txt.plain.slp": "$LADDER/aug_phoni/reference.txt.plain.slp",
  "augphoni/chr22-h10/reference.txt": "$LADDER/aug_phoni/reference.txt"
}}
EOF
$PY -m ms2exp.catalog hash --locations $T/toy-locations.json --out $T/toy-inputs.json
$PY -m ms2exp.catalog link --locations $T/toy-locations.json --root $T/inputs-root
$PY -m ms2exp.manifest --matrix toy --inputs $T/toy-inputs.json --ms2-repo $MS2_REPO --ms2-commit $COMMIT \
    --note "toy end-to-end run; mtb152 indexes built at the ms2 commit; reuses the 1d580b6 chr22-h10 min index" -o $T/manifest-toy.json
