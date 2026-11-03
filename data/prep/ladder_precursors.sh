#!/bin/bash
# RLBWT precursors of the A1 ladder references chr22-h1, h2, h5, h20 and h50 (precursors/chr22-hN/{bwt.heads,
# bwt.len,min}). These rungs were built by suffix sorting, not from an FMD: the text is every forward
# segment of reference.fa (data/prep/chr22_panel.py) and then their reverse complements, joined by % with a
# trailing newline; build-rlbwt (data/prep/build_rlbwt.cpp, a libdivsufsort64 suffix sort) gives text.bwt.heads
# and text.bwt.len, and TeraLCP -f rlbwt (BenLangmead/TeraTools af42783, branch teralcp-ms-minima) the minima
# text.min. Build build-rlbwt with
#   g++ -O3 -std=c++17 -I SDSL/include data/prep/build_rlbwt.cpp -L SDSL/lib -ldivsufsort64 -o build-rlbwt
# These are the build steps of ms-bakery bench/campaign-20261012/fig_rep_ladder_build.sbatch (Oct 27 2026);
# the index builds that followed them there are done by the freeze manifest instead.
# Usage: data/prep/ladder_precursors.sh REFERENCE_FA OUTDIR   (environment: BUILD_RLBWT, TERALCP, TMPDIR)
set -euo pipefail
REF=$1; OUT=$2
mkdir -p "$OUT"
cd "$OUT"
work=$(mktemp -d "${TMPDIR:-/tmp}/ladder-XXXXXX")
trap 'rm -rf "$work"' EXIT
python3 - "$REF" text.txt <<'PY'
import sys
seqs, cur = [], []
for line in open(sys.argv[1]):
    if line.startswith('>'):
        if cur:
            seqs.append(''.join(cur))
        cur = []
    else:
        cur.append(line.strip())
if cur:
    seqs.append(''.join(cur))
rc = lambda s: s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
with open(sys.argv[2], 'w') as f:
    f.write('%'.join(seqs + [rc(s) for s in seqs]) + '\n')
print('records', len(seqs), 'forward bases', sum(map(len, seqs)))
PY
$BUILD_RLBWT text.txt text
$TERALCP -f rlbwt -i text -t "$work/tmp" -ominima text.min
# text.txt is reproducible from the reference and is the largest file.
rm -f text.txt
sha256sum text.bwt.heads text.bwt.len text.min
