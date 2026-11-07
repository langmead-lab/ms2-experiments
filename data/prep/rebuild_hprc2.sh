#!/bin/bash
# Rebuilds the HPRC2 inputs that are not deposited, from public sources only, and checks each against
# the sha256 in data/inputs.tsv:
#   ropebwt3/human472/index.fmd                    Li's human472.fmr (Zenodo 14854401) converted to FMD
#   precursors/human472/{bwt.heads,bwt.len}        fmd2rlbwt on that FMD
#   precursors/human472/{min,phi}                  TeraLCP -f fmd -ominima -ophi on that FMD
#   ropebwt3/human472/index.fmd.{s8,s6}.ssa        ropebwt3 ssa -s8 and -s6 on that FMD
# The two inputs Li's record holds unchanged (index.fmd.len.gz, index.fmd.s11.ssa) come from
# python3 -m ms2exp.fetch, not from this script.
#
# Stages, run in this order (default: all of them):
#   tools   clone and build fmd2rlbwt (BenLangmead/runperm 2a3bf10) and TeraLCP (BenLangmead/TeraTools
#           tag ms2-precursors, commit 414f809, with lh3/ropebwt3 b408339 in its ropebwt3/ directory)
#   fmd     download human472.fmr.gz, check its md5 against Li's record, ropebwt3 build -i -do
#   runs    fmd2rlbwt
#   lcp     TeraLCP
#   ssa     ropebwt3 ssa -s8, then -s6
#   verify  python3 -m ms2exp.fetch --verify-only on the names this script writes
# Each stage checks its outputs' sha256 as soon as it finishes.
#
# Costs of the frozen run (Rockfish, 2026-10-01 and 02): fmd 3 min and 38 GiB after the download and gunzip;
# runs 4 min and 0.07 GiB; lcp about 13 h with 48 threads and 171 GiB of memory; each SSA about 10 h with 46 to 48
# threads. Disk: 8 GiB for the .fmr.gz, 14 GiB for the FMD, 147 GiB for the four precursors and 413 GiB
# for the two SSAs, plus TeraLCP's temporary files under WORK/tmp.
#
# Usage: data/prep/rebuild_hprc2.sh WORK [STAGE ...]
#   WORK holds the tool builds, the download and temporary files. Outputs go to DATA/<name>, where DATA
#   is $MS2EXP_DATA (as for ms2exp.fetch). Environment: MS2EXP_DATA, THREADS (default 48), ROPEBWT3 (default
#   the pinned binary $MS2EXP_DATA/tools/ropebwt3-cfad8475/ropebwt3 from python3 -m ms2exp.fetch
#   tools/ropebwt3-cfad8475/ropebwt3 tools/runtime/libgcc_s.so.1 tools/runtime/libgomp.so.1), FMD2RLBWT
#   and TERALCP (default the stage-tools builds under WORK), TERATOOLS_REPO and RUNPERM_REPO (where stage tools
#   clones from; default the GitHub repositories).
# Example, the cheap check of the FMD and the run files (about 30 min, mostly the download; 40 GB of memory):
#   data/prep/rebuild_hprc2.sh /tmp/hprc2 tools fmd runs
set -euo pipefail
WORK=$(mkdir -p "$1" && cd "$1" && pwd); shift
DATA=${MS2EXP_DATA:?set MS2EXP_DATA to the data directory}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
T=${THREADS:-48}
STAGES=${*:-tools fmd runs lcp ssa verify}
ROPEBWT3=${ROPEBWT3:-$DATA/tools/ropebwt3-cfad8475/ropebwt3}
FMD2RLBWT=${FMD2RLBWT:-$WORK/runperm/src/ms/tools/fmd2rlbwt/fmd2rlbwt}
TERALCP=${TERALCP:-$WORK/teratools/src/TeraLCP/TeraLCP}
export LD_LIBRARY_PATH=$DATA/tools/runtime${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
FMD=$DATA/ropebwt3/human472/index.fmd
PRE=$DATA/precursors/human472
LI=https://zenodo.org/records/14854401/files

check() {  # NAME: compare DATA/NAME with its sha256 in data/inputs.tsv
    local want got
    want=$(awk -F'\t' -v n="$1" '$1 == n {print $3}' "$REPO/data/inputs.tsv")
    got=$(sha256sum "$DATA/$1" | cut -d' ' -f1)
    if [ "$got" = "$want" ]; then echo "$1 equals the frozen file ($got)"
    else echo "$1 differs from the frozen file: $got, expected $want" >&2; return 1; fi
}

for stage in $STAGES; do
  echo "== $stage $(date -Is)"
  case $stage in
  tools)
    cd "$WORK"
    [ -d runperm ] || git clone -q "${RUNPERM_REPO:-https://github.com/BenLangmead/runperm.git}" runperm
    git -C runperm checkout -q 2a3bf102d22c12e8c7d4dc9b3bf2f53818c1629d
    make -C runperm/src/ms/tools/fmd2rlbwt
    [ -d teratools ] || git clone -q "${TERATOOLS_REPO:-https://github.com/BenLangmead/TeraTools.git}" teratools
    git -C teratools checkout -q ms2-precursors
    test "$(git -C teratools rev-parse HEAD)" = 414f809dfd34c79547a28ff2a4c58e4e505ce94f
    # TeraTools' Makefiles clone lh3/ropebwt3 at its current head when ropebwt3/ is missing; pin it first.
    [ -d teratools/ropebwt3 ] || git clone -q https://github.com/lh3/ropebwt3 teratools/ropebwt3
    git -C teratools/ropebwt3 checkout -q b408339d0fb08fe0f4da5140507a6b8b854fa12f
    make -C teratools/ropebwt3
    make -C teratools/src thirdparty
    make -C teratools/src/TeraLCP TeraLCP
    sha256sum "$FMD2RLBWT" "$TERALCP"
    ;;
  fmd)
    mkdir -p "$(dirname "$FMD")"; cd "$WORK"
    curl -fsSL -C - --retry 5 -o human472.fmr.gz "$LI/human472.fmr.gz?download=1"
    echo "5304b0bce6ae545ffd4aee7fefdb6ea0  human472.fmr.gz" | md5sum -c -
    gzip -dc human472.fmr.gz > human472.fmr
    "$ROPEBWT3" build -t 8 -i human472.fmr -do "$FMD.partial"
    mv "$FMD.partial" "$FMD"; rm -f human472.fmr
    check ropebwt3/human472/index.fmd
    ;;
  runs)
    mkdir -p "$PRE"
    "$FMD2RLBWT" "$FMD" "$PRE/fmd" > "$WORK/fmd2rlbwt.txt"
    mv "$PRE/fmd.bwt.heads" "$PRE/bwt.heads"; mv "$PRE/fmd.bwt.len" "$PRE/bwt.len"
    check precursors/human472/bwt.heads; check precursors/human472/bwt.len
    ;;
  lcp)
    # TeraLCP 414f809 numbers the minima by the FMD's maximal runs, the numbering of bwt.heads and
    # bwt.len, so they need no renumbering (fmd2rlbwt --minima-in is only for older TeraLCP builds).
    mkdir -p "$PRE" "$WORK/tmp"
    /usr/bin/time -v "$TERALCP" -f fmd -i "$FMD" -t "$WORK/tmp/t" -p "$T" \
        -ominima "$PRE/min.partial" -ophi "$PRE/phi.partial" -v time > "$WORK/teralcp.stdout" 2> "$WORK/teralcp.stderr"
    mv "$PRE/min.partial" "$PRE/min"; mv "$PRE/phi.partial" "$PRE/phi"
    check precursors/human472/min; check precursors/human472/phi
    ;;
  ssa)
    for s in 8 6; do
      "$ROPEBWT3" ssa -t "$T" -s "$s" -o "$FMD.s$s.ssa.partial" "$FMD"
      mv "$FMD.s$s.ssa.partial" "$FMD.s$s.ssa"
      check ropebwt3/human472/index.fmd.s$s.ssa
    done
    ;;
  verify)
    cd "$REPO"
    python3 -m ms2exp.fetch --data "$DATA" --verify-only ropebwt3/human472/index.fmd \
        precursors/human472/bwt.heads precursors/human472/bwt.len precursors/human472/min \
        precursors/human472/phi ropebwt3/human472/index.fmd.s8.ssa ropebwt3/human472/index.fmd.s6.ssa
    ;;
  *) echo "unknown stage $stage (tools, fmd, runs, lcp, ssa or verify)" >&2; exit 1 ;;
  esac
done
