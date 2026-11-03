#!/usr/bin/env python3
"""Convert a (gzipped) FASTA assembly into ACGT-only records.

Rule: uppercase every base, split each input record at every run of
characters other than A, C, G, T, and drop pieces shorter than 1 base.
Each piece is named <origname>_<start>_<end>, where origname is the first
whitespace-delimited word of the input header and start/end are 0-based,
half-open coordinates of the piece within the input record. Output is a
two-line FASTA (header, then the whole sequence on one line).

Also writes offsets.u64: little-endian uint64 cumulative start offsets of the
output records (n_records + 1 values, first 0, last the total output bases),
and a JSON file of statistics.

Memory use is bounded by the longest output piece (held in a bytearray).

Usage: python3 -I acgt_split.py IN.fa[.gz] OUT.fa OUT.offsets.u64 OUT.stats.json
"""

import gzip
import hashlib
import json
import re
import struct
import sys

ACGT_RUN = re.compile(rb"[ACGT]+")


def main():
    inp, out_fa, out_off, out_stats = sys.argv[1:5]
    opener = gzip.open if inp.endswith(".gz") else open
    st = {
        "input_records": 0,
        "input_bases": 0,
        "input_acgt_bases": 0,
        "input_lowercase_bases": 0,
        "non_acgt_bases_dropped": 0,
        "non_acgt_by_char": {},
        "output_records": 0,
        "output_bases": 0,
        "input_records_with_non_acgt": 0,
        "input_records_fully_dropped": 0,
        "longest_output_record": 0,
        "shortest_output_record": None,
    }
    sha = hashlib.sha256()
    fout = open(out_fa, "wb")
    foff = open(out_off, "wb")
    foff.write(struct.pack("<Q", 0))

    name = None
    pos = 0            # position within current input record
    piece = bytearray()
    piece_start = 0
    rec_pieces = 0
    rec_had_non = False

    def flush():
        nonlocal piece, rec_pieces
        if len(piece) >= 1:
            hdr = b">%s_%d_%d\n" % (name, piece_start, piece_start + len(piece))
            fout.write(hdr); sha.update(hdr)
            fout.write(piece); sha.update(piece)
            fout.write(b"\n"); sha.update(b"\n")
            n = len(piece)
            st["output_records"] += 1
            st["output_bases"] += n
            foff.write(struct.pack("<Q", st["output_bases"]))
            st["longest_output_record"] = max(st["longest_output_record"], n)
            s = st["shortest_output_record"]
            st["shortest_output_record"] = n if s is None else min(s, n)
            rec_pieces += 1
        piece = bytearray()

    def end_record():
        if name is None:
            return
        flush()
        if rec_had_non:
            st["input_records_with_non_acgt"] += 1
        if rec_pieces == 0:
            st["input_records_fully_dropped"] += 1

    nonacgt = st["non_acgt_by_char"]
    with opener(inp, "rb") as f:
        for line in f:
            line = line.rstrip(b"\r\n")
            if line.startswith(b">"):
                end_record()
                name = line[1:].split()[0]
                st["input_records"] += 1
                pos = 0
                piece = bytearray()
                piece_start = 0
                rec_pieces = 0
                rec_had_non = False
                continue
            if not line:
                continue
            L = len(line)
            st["input_bases"] += L
            up = line.upper()
            if up != line:
                st["input_lowercase_bases"] += sum(1 for c in line if 97 <= c <= 122)
            cur = 0
            for m in ACGT_RUN.finditer(up):
                a, b = m.span()
                if a > cur:
                    # non-ACGT run at [cur, a) ends the current piece
                    for c in up[cur:a]:
                        k = chr(c)
                        nonacgt[k] = nonacgt.get(k, 0) + 1
                    st["non_acgt_bases_dropped"] += a - cur
                    rec_had_non = True
                    flush()
                if not piece:
                    piece_start = pos + a
                piece += up[a:b]
                st["input_acgt_bases"] += b - a
                cur = b
            if cur < L:
                for c in up[cur:L]:
                    k = chr(c)
                    nonacgt[k] = nonacgt.get(k, 0) + 1
                st["non_acgt_bases_dropped"] += L - cur
                rec_had_non = True
                flush()
            pos += L
    end_record()
    fout.close()
    foff.close()
    st["reads_fa_sha256"] = sha.hexdigest()
    assert st["output_bases"] == st["input_acgt_bases"]
    assert st["input_bases"] == st["input_acgt_bases"] + st["non_acgt_bases_dropped"]
    with open(out_stats, "w") as g:
        json.dump(st, g, indent=1, sort_keys=True)
        g.write("\n")
    print(json.dumps(st, sort_keys=True))


if __name__ == "__main__":
    main()
