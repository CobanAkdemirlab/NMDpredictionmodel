#!/usr/bin/env python3
import argparse, re, os
import pandas as pd
import numpy as np
import pyranges as pr
from pyfaidx import Fasta
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# ------------------ helpers ------------------

def harmonize_chr(df, col="Chromosome", want_chr_prefix=True):
    s = df[col].astype(str)
    has_chr = s.str.startswith("chr").any()
    if want_chr_prefix and not has_chr:
        df[col] = "chr" + s.str.replace("^chr", "", regex=True)
    if not want_chr_prefix and has_chr:
        df[col] = s.str.replace("^chr", "", regex=True)
    return df

def strip_version(x):
    if pd.isna(x): return x
    return str(x).split(".", 1)[0]

def split_txnames(val):
    if pd.isna(val): return []
    parts = re.split(r"[;,|\s]+", str(val).strip())
    return [p for p in parts if p]

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

def load_optimal_codons(path: str) -> set:
    s = set()
    with open(path) as f:
        for line in f:
            c = line.strip().upper().replace("U", "T")
            if len(c) == 3 and all(b in "ACGT" for b in c):
                s.add(c)
    if not s:
        raise SystemExit("ERROR: optimal codon set is empty. Check --optimal_codons.")
    return s

# ------------------ load variants ------------------

def load_variants_csv(path):
    """
    Required columns: contig, position, txnames.
    Optional: variantID -> PTC_ID
    Returns v_orig (PTC-level), v_exp (exploded per transcript)
    """
    v = pd.read_csv(path)
    required = {"contig", "position", "txnames"}
    missing = required - set(v.columns)
    if missing:
        raise SystemExit(f"ERROR: Missing required columns: {missing}")

    v = v.rename(columns={"contig": "contig", "position": "Position"})
    if "variantID" in v.columns and "PTC_ID" not in v.columns:
        v = v.rename(columns={"variantID": "PTC_ID"})
    if "PTC_ID" not in v.columns:
        v["PTC_ID"] = v.index.astype(str)

    v["Chromosome"] = v["contig"].astype(str)
    v["Position"]   = v["Position"].astype(int)
    v = harmonize_chr(v, "Chromosome", want_chr_prefix=True)

    rows = []
    for _, r in v.iterrows():
        txs = [strip_version(t) for t in split_txnames(r["txnames"])]
        if not txs:
            rows.append([r["PTC_ID"], r["Chromosome"], r["Position"], np.nan])
        else:
            for t in txs:
                rows.append([r["PTC_ID"], r["Chromosome"], r["Position"], t])
    v_exp = pd.DataFrame(rows, columns=["PTC_ID","Chromosome","Position","TxName"])
    return v[["PTC_ID","Chromosome","Position"]].drop_duplicates(), v_exp

# ------------------ GTF -> CDS model per transcript ------------------

def load_cds_from_gtf(gtf_path):
    """
    Returns cds_by_tx: dict tx -> DataFrame(Chromosome, Start, End, Strand), where
    Start/End are 0-based half-open (PyRanges convention).
    """
    gtf = pr.read_gtf(gtf_path)
    gdf = gtf.df.copy()
    if "transcript_id" not in gdf.columns:
        raise SystemExit("ERROR: GTF missing transcript_id attribute.")

    gdf["TxName"] = gdf["transcript_id"].astype(str).map(strip_version)
    cds = gdf[gdf.Feature == "CDS"][["Chromosome","Start","End","Strand","TxName"]].copy()
    cds = harmonize_chr(cds, "Chromosome", want_chr_prefix=True)

    cds_by_tx = {}
    for tx, d in cds.groupby("TxName", sort=False):
        d = d.sort_values(["Chromosome","Start","End"]).reset_index(drop=True)
        cds_by_tx[tx] = d
    return cds_by_tx

# ------------------ codon optimality ------------------

def cds_segments_in_coding_order(df_cds_tx: pd.DataFrame):
    strand = df_cds_tx["Strand"].iloc[0]
    segs = [(r["Chromosome"], int(r["Start"]), int(r["End"])) for _, r in df_cds_tx.iterrows()]
    return segs if strand == "+" else list(reversed(segs))

def genomic_pos_to_cds_offset_nt(pos_1based: int, df_cds_tx: pd.DataFrame):
    """
    Map genomic POS (1-based, VCF-style) -> CDS offset (0-based nt in coding orientation).
    Returns None if POS not inside CDS.
    """
    pos0 = pos_1based - 1
    strand = df_cds_tx["Strand"].iloc[0]
    offset = 0
    for chrom, s0, e0 in cds_segments_in_coding_order(df_cds_tx):
        if s0 <= pos0 < e0:
            if strand == "+":
                return offset + (pos0 - s0)
            else:
                return offset + ((e0 - 1) - pos0)
        offset += (e0 - s0)
    return None

def fetch_spliced_cds_seq(genome: Fasta, df_cds_tx: pd.DataFrame):
    strand = df_cds_tx["Strand"].iloc[0]
    chunks = []
    for chrom, s0, e0 in cds_segments_in_coding_order(df_cds_tx):
        chunks.append(genome[chrom][s0:e0].seq.upper())
    seq = "".join(chunks)
    if strand == "-":
        seq = revcomp(seq)
    seq = seq[: (len(seq)//3)*3]  # full codons only
    return seq

def codons_from_seq(seq: str):
    seq = seq.upper().replace("U","T")
    return [seq[i:i+3] for i in range(0, len(seq)-2, 3) if "N" not in seq[i:i+3]]

def optimal_fraction(codons, optimal_codons):
    if not codons:
        return np.nan
    return sum(c in optimal_codons for c in codons) / len(codons)

def optimal_fraction_window_pm_nt(codons, ptc_offset_nt, optimal_codons, window_nt=100):
    if not codons or ptc_offset_nt is None:
        return np.nan
    cds_len_nt = len(codons) * 3
    lo_nt = max(0, ptc_offset_nt - window_nt)
    hi_nt = min(cds_len_nt - 1, ptc_offset_nt + window_nt)
    lo_c = lo_nt // 3
    hi_c = hi_nt // 3
    return optimal_fraction(codons[lo_c:hi_c+1], optimal_codons)

# ------------------ parallel workers ------------------

def _work_chunk(args):
    chunk_df, cds_by_tx, fasta_path, optimal_codons, window_nt = args
    genome = Fasta(fasta_path)

    tx_cache = {}  # tx -> codon list
    out_rows = []

    for _, r in chunk_df.iterrows():
        pid = r["PTC_ID"]
        chrom = r["Chromosome"]
        pos = int(r["Position"])
        tx = r["TxName"]

        if pd.isna(tx) or tx not in cds_by_tx:
            out_rows.append([pid, chrom, pos, tx, np.nan, np.nan, "no_tx_or_no_cds"])
            continue

        df_cds_tx = cds_by_tx[tx]

        # optional sanity: transcript has CDS on this chromosome
        if not (df_cds_tx["Chromosome"] == chrom).any():
            out_rows.append([pid, chrom, pos, tx, np.nan, np.nan, "chr_mismatch"])
            continue

        if tx not in tx_cache:
            cds_seq = fetch_spliced_cds_seq(genome, df_cds_tx)
            tx_cache[tx] = codons_from_seq(cds_seq)

        codons = tx_cache[tx]
        whole = optimal_fraction(codons, optimal_codons)

        ptc_off = genomic_pos_to_cds_offset_nt(pos, df_cds_tx)
        win = optimal_fraction_window_pm_nt(codons, ptc_off, optimal_codons, window_nt=window_nt)

        status = "ok" if ptc_off is not None else "ptc_not_in_cds"
        out_rows.append([pid, chrom, pos, tx, whole, win, status])

    return pd.DataFrame(out_rows, columns=[
        "PTC_ID","Chromosome","Position","TxName",
        "CodonOptimalityFraction_CDS",
        "CodonOptimalityFraction_PTCpm100nt",
        "status"
    ])

def parallel_compute(v_exp, cds_by_tx, fasta_path, optimal_codons, window_nt=100, threads=4, chunksize=50_000):
    if v_exp.empty:
        return pd.DataFrame(columns=[
            "PTC_ID","Chromosome","Position","TxName",
            "CodonOptimalityFraction_CDS","CodonOptimalityFraction_PTCpm100nt","status"
        ])

    chunks = [v_exp.iloc[i:i+chunksize].copy() for i in range(0, len(v_exp), chunksize)]
    tasks = [(c, cds_by_tx, fasta_path, optimal_codons, window_nt) for c in chunks]

    out = []
    with ProcessPoolExecutor(max_workers=max(1, threads)) as ex:
        futures = [ex.submit(_work_chunk, t) for t in tasks]
        for i, f in enumerate(as_completed(futures), 1):
            out.append(f.result())
            if i % 5 == 0:
                print(f"    … completed {i}/{len(tasks)} chunks")
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()

# ------------------ main ------------------

def main():
    ap = argparse.ArgumentParser(
        description="Compute codon optimality fraction per transcript (txnames-driven): whole CDS + ±100nt window around PTC."
    )
    ap.add_argument("--variants", required=True, help="CSV with contig, position, txnames, [variantID]")
    ap.add_argument("--gtf", required=True, help="Gencode GTF (gz OK)")
    ap.add_argument("--fasta", default="/Users/jschmidt3/Iman_visualizations/efficient_motif/hg38.fa",
                    help="Genome FASTA matching the GTF (default: your hg38.fa)")
    ap.add_argument("--optimal_codons", required=True, help="Text file: one optimal codon per line (e.g., GCC)")
    ap.add_argument("--out", required=True, help="Output TSV[.gz]")
    ap.add_argument("--window_nt", type=int, default=100, help="Half-window in nt around PTC (default 100)")
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 4, help="Workers (default all cores)")
    ap.add_argument("--aggregate", choices=["none","max","mean","min"], default="none",
                    help="Aggregate across multiple txnames per PTC (default none)")
    args = ap.parse_args()

    print("[*] Loading variants…")
    v_orig, v_exp = load_variants_csv(args.variants)
    print(f"    - PTC rows: {len(v_orig):,}")
    print(f"    - Exploded rows (PTC x txnames): {len(v_exp):,}")

    print("[*] Loading optimal codons…")
    optimal_codons = load_optimal_codons(args.optimal_codons)
    print(f"    - Optimal codons loaded: {len(optimal_codons)}")

    print("[*] Loading CDS from GTF…")
    cds_by_tx = load_cds_from_gtf(args.gtf)
    print(f"    - Transcripts with CDS: {len(cds_by_tx):,}")

    print(f"[*] Computing codon optimality (threads={args.threads})…")
    res = parallel_compute(
        v_exp=v_exp,
        cds_by_tx=cds_by_tx,
        fasta_path=args.fasta,
        optimal_codons=optimal_codons,
        window_nt=args.window_nt,
        threads=args.threads
    )

    # Optional aggregation across txnames per PTC
    if args.aggregate != "none" and len(res):
        aggfunc = {"max": "max", "mean": "mean", "min": "min"}[args.aggregate]
        res = (res.groupby("PTC_ID", as_index=False)
                 .agg({
                     "CodonOptimalityFraction_CDS": aggfunc,
                     "CodonOptimalityFraction_PTCpm100nt": aggfunc
                 }))

    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    if outp.suffix.endswith("gz"):
        res.to_csv(outp, sep="\t", index=False, compression="gzip")
    else:
        res.to_csv(outp, sep="\t", index=False)

    print(f"[*] Done. Wrote {len(res):,} rows to {outp}")
    print("[*] Status counts:")
    print(res["status"].value_counts(dropna=False).to_string())

if __name__ == "__main__":
    main()
