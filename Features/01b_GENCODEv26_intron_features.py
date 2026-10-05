#!/usr/bin/env python3
"""
01b_GENCODEv26_intron_features.py

Generate transcript-level CDS and UTR intron features from
GENCODE v26 exon and CDS annotations.

Features:
    cdsseq.introns
    threeUTR.introns
    fiveUTR.introns

Definitions:
    cdsseq.introns
        TRUE when the CDS, including the stop codon, spans >1 exon.

    threeUTR.introns
        TRUE when the annotated 3'UTR spans >1 exon.

    fiveUTR.introns
        TRUE when the annotated 5'UTR spans >1 exon.

Transcript strand is used to distinguish the 5' and 3' sides
of the CDS.

Reference:
    GENCODE v26 primary assembly
    GRCh38 / hg38

"""
import argparse, gzip, re
from collections import Counter

FLAGS = {"cdsseq.introns": '"There is a cdsseq intron"',
         "threeUTR.introns": '"There is a 3UTR intron"',
         "fiveUTR.introns": '"There is a 5UTR intron"'}


def split_raw(line):
    out, cur, q = [], [], False
    for ch in line:
        if ch == '"':
            q = not q
        if ch == ',' and not q:
            out.append(''.join(cur)); cur = []
        else:
            cur.append(ch)
    out.append(''.join(cur))
    return out


def unq(x):
    return x[1:-1] if len(x) >= 2 and x[0] == x[-1] == '"' else x


ap = argparse.ArgumentParser()
ap.add_argument("--table", required=True)
ap.add_argument("--gtf", required=True)
ap.add_argument("--output", required=True)
a = ap.parse_args()

text = open(a.table, newline="").read()
eol = "\r\n" if "\r\n" in text[:100000] else "\n"
lines = text.split(eol)
header = split_raw(lines[0]); names = [unq(h) for h in header]
tx_i = names.index("txnames")
idx = {c: names.index(c) for c in FLAGS if c in names}
want = {unq(split_raw(l)[tx_i]) for l in lines[1:] if l}

exons, cds, strand = {}, {}, {}
op = gzip.open if a.gtf.endswith(".gz") else open
with op(a.gtf, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 9 or p[2] not in ("exon", "CDS", "stop_codon"):
            continue
        m = re.search(r'transcript_id "([^"]+)"', p[8])
        if not m or m.group(1) not in want:
            continue
        tx = m.group(1); s, e = int(p[3]), int(p[4]); strand[tx] = p[6]
        (exons if p[2] == "exon" else cds).setdefault(tx, []).append((s, e))


def flags(tx):
    ex, c = exons.get(tx), cds.get(tx)
    if not ex or not c:
        return None
    lo, hi = min(s for s, _ in c), max(e for _, e in c)          # CDS incl. stop codon
    n_cds = sum(1 for s, e in ex if s <= hi and e >= lo)
    n_up = sum(1 for s, e in ex if s < lo)                        # exons with sequence left of CDS
    n_dn = sum(1 for s, e in ex if e > hi)                        # exons with sequence right of CDS
    n5, n3 = (n_up, n_dn) if strand[tx] == "+" else (n_dn, n_up)
    return {"cdsseq.introns": n_cds > 1, "threeUTR.introns": n3 > 1, "fiveUTR.introns": n5 > 1}


out = [lines[0]]; changed = Counter(); missing = 0; yes = Counter()
for l in lines[1:]:
    if not l:
        out.append(l); continue
    fl = split_raw(l)
    fx = flags(unq(fl[tx_i]))
    if fx is None:
        missing += 1
    for c, i in idx.items():
        new = "NA" if fx is None or not fx[c] else FLAGS[c]
        if new != "NA":
            yes[c] += 1
        if (fl[i] if fl[i] else "NA") != new:
            changed[c] += 1
        fl[i] = new
    out.append(",".join(fl))

open(a.output, "w", newline="").write(eol.join(out))
n = sum(1 for l in lines[1:] if l)
print(f"variants: {n} | transcripts not found in GTF: {missing}")
for c in idx:
    print(f"  {c:18s} flagged {yes[c]:5d} | changed {changed[c]:5d}")
print(f"Wrote {a.output}")
