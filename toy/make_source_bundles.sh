#!/bin/bash
# Bundle the ms2 source of a freeze commit, with its Orbit and TeraTools submodules, so that build
# boxes build it without access to the private langmead-lab/ms2 repository.
# Usage: make_source_bundles.sh MS2_CHECKOUT MS2_COMMIT OUT_DIR [HASH_DIR]
# MS2_CHECKOUT must have MS2_COMMIT checked out, clean, with both submodules at the commits it
# records (git submodule update --init). Writes OUT_DIR/{ms2,orbit,teratools}.bundle, each holding
# its commit as HEAD with its history, and OUT_DIR/source-locations.json, which maps the logical
# names source/ms2-<c7>/{ms2,orbit,teratools}.bundle to them for catalog.py hash and link. HASH_DIR
# (default OUT_DIR) is where the bundles will be when catalog.py runs, as when OUT_DIR is made on the
# Mac and copied to Rockfish. The build-tools job of an all-AWS manifest (ms2exp.manifest
# --ms2-source bundle, the default) stages the three bundles as inputs.
set -euo pipefail
C=$(cd "$1" && pwd)
COMMIT=$2
mkdir -p "$3"
OUT=$(cd "$3" && pwd)
H=${4:-$OUT}
[ ${#COMMIT} -eq 40 ] || { echo "MS2_COMMIT must be a full 40-hex commit" >&2; exit 1; }
[ "$(git -C "$C" rev-parse HEAD)" = "$COMMIT" ] || { echo "$C is not at $COMMIT" >&2; exit 1; }
git -C "$C" diff --quiet HEAD || { echo "$C has uncommitted changes" >&2; exit 1; }
git -C "$C" bundle create -q "$OUT/ms2.bundle" HEAD
for m in orbit teratools; do
    want=$(git -C "$C" ls-tree HEAD "$m" | awk '{print $3}')
    have=$(git -C "$C/$m" rev-parse HEAD)
    [ "$want" = "$have" ] || { echo "$m is at $have; $COMMIT records $want" >&2; exit 1; }
    git -C "$C/$m" bundle create -q "$OUT/$m.bundle" HEAD
done
for b in ms2 orbit teratools; do
    git bundle list-heads "$OUT/$b.bundle"
done
S=source/ms2-${COMMIT:0:7}
cat > "$OUT/source-locations.json" <<EOF
{"locations": {
  "$S/ms2.bundle": "$H/ms2.bundle",
  "$S/orbit.bundle": "$H/orbit.bundle",
  "$S/teratools.bundle": "$H/teratools.bundle"
}}
EOF
