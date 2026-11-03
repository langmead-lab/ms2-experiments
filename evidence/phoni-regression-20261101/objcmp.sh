#!/bin/bash
# Compares a freeze binary with a rebuild: file sha256, section sizes, disassembly and read-only strings.
# Usage: objcmp.sh FREEZE REBUILD
a=$1; b=$2
echo "freeze:  $(sha256sum < $a | cut -c1-64)  $a"
echo "rebuild: $(sha256sum < $b | cut -c1-64)  $b"
echo "section sizes (name freeze rebuild):"
join <(readelf -SW $a | awk '$1 ~ /^\[/ && $2 ~ /^\./ {print $2, $6} $1 ~ /^\[/ && $3 ~ /^\./ {print $3, $7}' | sort) \
     <(readelf -SW $b | awk '$1 ~ /^\[/ && $2 ~ /^\./ {print $2, $6} $1 ~ /^\[/ && $3 ~ /^\./ {print $3, $7}' | sort) |
  awk '{printf "  %-22s %s %s%s\n", $1, $2, $3, ($2 == $3 ? "" : "  differs")}'
dis() { objdump -d --no-show-raw-insn "$1" | grep -v 'file format' | sed -E 's/^ *[0-9a-f]+:\t//'; }
echo "disassembly lines: $(dis $a | wc -l) and $(dis $b | wc -l)"
echo "differing lines with only instruction addresses removed: $(diff <(dis $a) <(dis $b) | grep -c '^[<>]')"
echo "differing lines with every hex number of 4 or more digits masked: $(diff <(dis $a | sed -E 's/[0-9a-f]{4,}/H/g') <(dis $b | sed -E 's/[0-9a-f]{4,}/H/g') | grep -c '^[<>]')"
echo "strings found in only one binary (strings -a -n 8):"
diff <(strings -a -n 8 $a | sort -u) <(strings -a -n 8 $b | sort -u) | grep '^[<>]' | sed 's/^/  /'
