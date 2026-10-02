#!/usr/bin/env python3
"""
13_isoform_count_features.py

Purpose
-------
Add the `isoform_count` feature to a variant feature table.

isoform_count = number of distinct GENCODE v26 transcripts annotated for the
variant's gene (all transcript records of the gene, regardless of biotype).
This is the definition used for the TOPMed training table, so TOPMed, gnomAD,
ClinVar and GREGoR must all be built with this script and the same GTF.

Method
------
1. Parse the GTF `transcript` records; collect the unique transcript_ids of each
   gene_id (gene version suffix stripped, e.g. ENSG00000123456.7 -> ENSG00000123456).
2. Map the table's gene ID column (default `ensembl_gene_id`, version suffix
   stripped on the table side too) to that count.
3. Genes that cannot be mapped receive 1, exactly as in the TOPMed training table.
   The number and share of such variants is reported, and the run stops if the
   share exceeds --max_missing_frac (default 5%), which usually means the wrong GTF
   or a wrong gene ID column.

Input / output
--------------
The table is read in full and written to --output with an added (or replaced)
`isoform_count` column. The input file is never overwritten unless you pass
--overwrite.

Usage
-----
python 13_isoform_count_features.py \
    --table  clinvar_base.csv \
    --gtf    gencode.v26.primary_assembly.annotation.gtf.gz \
    --output clinvar_base_isoform.csv

Notes
-----
* Replaces the earlier one-off `add_isoform_counts.py`, which had hard-coded
  paths, wrote its result over the input file, did not strip version suffixes
  from the table's gene IDs (versioned IDs would all have silently become 1),
  and did not report how many variants were filled.
"""
import argparse
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

GENE_RE = re.compile(r'gene_id "([^"]+)"')
TX_RE = re.compile(r'transcript_id "([^"]+)"')


def strip_version(x):
    return x if pd.isna(x) else str(x).split(".", 1)[0]


def count_isoforms(gtf_path: str) -> dict:
    """gene_id (no version) -> number of distinct transcript_ids."""
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    tx_by_gene = defaultdict(set)
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            g, t = GENE_RE.search(fields[8]), TX_RE.search(fields[8])
            if g and t:
                tx_by_gene[g.group(1).split(".", 1)[0]].add(t.group(1))
    return {g: len(t) for g, t in tx_by_gene.items()}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, type=Path, help="variant feature table (CSV)")
    ap.add_argument("--gtf", required=True, help="GENCODE v26 GTF (gz OK)")
    ap.add_argument("--output", required=True, type=Path, help="output CSV (table + isoform_count)")
    ap.add_argument("--id_col", default="ensembl_gene_id", help="gene ID column in the table (default ensembl_gene_id)")
    ap.add_argument("--max_missing_frac", type=float, default=0.05,
                    help="stop if more than this share of variants cannot be mapped (default 0.05)")
    ap.add_argument("--overwrite", action="store_true", help="allow --output to equal --table")
    args = ap.parse_args()

    if args.output.resolve() == args.table.resolve() and not args.overwrite:
        sys.exit("ERROR: --output is the same file as --table. Choose another path or pass --overwrite.")

    print("[*] Parsing GTF for isoform counts...")
    counts = count_isoforms(args.gtf)
    print(f"    genes with transcripts: {len(counts):,}")

    print(f"[*] Loading {args.table}")
    df = pd.read_csv(args.table, low_memory=False)
    print(f"    shape: {df.shape}")
    if args.id_col not in df.columns:
        sys.exit(f"ERROR: column '{args.id_col}' not found. Columns containing 'gene': "
                 f"{[c for c in df.columns if 'gene' in c.lower()][:10]}")

    gid = df[args.id_col].map(strip_version)
    mapped = gid.map(counts)
    n_missing = int(mapped.isna().sum())
    frac = n_missing / len(df)
    print(f"[*] variants whose gene could not be mapped (filled with 1): {n_missing:,} ({frac:.2%})")
    if n_missing:
        print(f"    examples: {sorted(set(gid[mapped.isna()].dropna().astype(str)))[:5]}")
    if frac > args.max_missing_frac:
        sys.exit(f"ERROR: {frac:.1%} unmapped exceeds --max_missing_frac {args.max_missing_frac:.0%}: "
                 "check the GTF version and the gene ID column.")

    df["isoform_count"] = mapped.fillna(1).astype(int)
    print(f"    isoform_count: min {df['isoform_count'].min()} | median {df['isoform_count'].median():.0f} "
          f"| max {df['isoform_count'].max()} | ==1: {(df['isoform_count'] == 1).sum():,}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"[*] Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
