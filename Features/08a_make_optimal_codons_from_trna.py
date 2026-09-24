#!/usr/bin/env python3
"""
08a_make_optimal_codons_from_trna.py

Purpose
-------
Generate the optimal-codon reference used for codon-optimality feature
calculation.

Optimal codons were defined following the approach used by
Lindeboom et al. (2016): for each amino acid, the codon(s) corresponding
to the highest number of human tRNA genes were classified as optimal.

Reference resource
------------------
Human tRNA annotation:
    UCSC Table Browser
    Genome assembly: hg38
    Track: GtRNAdb-based tRNA genes

Repository input:
    reference/codon_optimality/hg38_tRNAs.tsv

The input table must contain:
    aa  - amino acid
    ac  - anticodon

Method
------
1. Reverse-complement each tRNA anticodon to obtain the corresponding
   DNA codon.
2. Count tRNA gene copies for each amino-acid/codon combination.
3. For each amino acid, retain the codon(s) with the maximum tRNA gene
   copy number.

The manuscript analysis used:
    --mode per_aa_max

Output
------
Plain-text file containing one optimal codon per line:

    reference/codon_optimality/optimal_codons.txt

This output is subsequently used by:

    08b_codon_optimality_features.py

References
----------
Lindeboom RGH, Supek F, Lehner B.
The rules and impact of nonsense-mediated mRNA decay in human cancers.
Nature Genetics. 2016.

Supek F, Minana B, Valcarcel J, Gabaldon T, Lehner B.
Synonymous mutations frequently act as driver mutations in human cancers.
Cell. 2014.
"""

import argparse
import pandas as pd
from pathlib import Path

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

def main():
    ap = argparse.ArgumentParser(description="Make optimal codon list from hg38 tRNA gene copy numbers.")
    ap.add_argument("--trna_tsv", required=True, help="UCSC tRNAs table export TSV (must include aa and ac columns)")
    ap.add_argument("--aa_col", default="aa", help="Column name for amino acid (default: aa)")
    ap.add_argument("--ac_col", default="ac", help="Column name for anticodon (default: ac)")
    ap.add_argument("--mode", choices=["per_aa_max","top_quantile"], default="per_aa_max",
                    help="per_aa_max: pick codon(s) with max tRNA copies per amino acid; "
                         "top_quantile: pick codons in top quantile by tRNA copies (global).")
    ap.add_argument("--quantile", type=float, default=0.80, help="Quantile for top_quantile mode (default 0.80)")
    ap.add_argument("--out", required=True, help="Output optimal_codons.txt")
    args = ap.parse_args()

    df = pd.read_csv(args.trna_tsv, sep="\t")
    for col in (args.aa_col, args.ac_col):
        if col not in df.columns:
            raise SystemExit(f"Missing column '{col}' in {args.trna_tsv}. Columns: {list(df.columns)}")

    df = df[[args.aa_col, args.ac_col]].dropna().copy()
    df[args.ac_col] = df[args.ac_col].astype(str).str.upper().str.replace("U","T")
    df[args.aa_col] = df[args.aa_col].astype(str)

    # anticodon -> codon (DNA): codon is reverse-complement of anticodon
    df["codon"] = df[args.ac_col].apply(lambda x: revcomp(x) if len(x) == 3 else None)
    df = df[df["codon"].notna()].copy()

    # count tRNA gene copies by (aa, codon)
    counts = (df.groupby([args.aa_col, "codon"])
                .size().reset_index(name="trna_gene_copies"))

    if args.mode == "per_aa_max":
        # for each amino acid, pick codon(s) with maximum copy number
        counts["max_for_aa"] = counts.groupby(args.aa_col)["trna_gene_copies"].transform("max")
        optimal = counts[counts["trna_gene_copies"] == counts["max_for_aa"]].copy()
    else:
        # global top quantile by copy number
        thr = counts["trna_gene_copies"].quantile(args.quantile)
        optimal = counts[counts["trna_gene_copies"] >= thr].copy()

    codons = sorted(set(optimal["codon"].tolist()))
    with open(args.out, "w") as f:
        for c in codons:
            f.write(c + "\n")

    print(f"[✓] Wrote {len(codons)} optimal codons to {args.out}")
    print(optimal.sort_values("trna_gene_copies", ascending=False).head(10).to_string(index=False))

if __name__ == "__main__":
    main()
