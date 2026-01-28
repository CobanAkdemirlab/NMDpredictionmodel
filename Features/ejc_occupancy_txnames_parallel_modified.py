#!/usr/bin/env python3
import argparse, re, os
import pandas as pd
import numpy as np
import pyranges as pr
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
    # ENST00000335137.4 -> ENST00000335137
    if pd.isna(x): return x
    return str(x).split('.', 1)[0]

def split_txnames(val):
    """Split multi-valued txnames like 'ENST...;ENST...' or comma/space/pipe separated."""
    if pd.isna(val): return []
    parts = re.split(r"[;,|\s]+", str(val).strip())
    parts = [p for p in parts if p]
    return parts

# ------------------ load inputs ------------------

def load_variants_csv(path):
    """
    Always read variants as CSV. Required columns: contig, position, txnames.
    Optional: variantID (used as PTC_ID).
    Returns:
      v_orig: original per-PTC rows (with PTC_ID assigned)
      v_exp: exploded rows per (PTC_ID, TxName)
    """
    v = pd.read_csv(path)

    required = {"contig", "position", "txnames"}
    missing = required - set(v.columns)
    if missing:
        raise SystemExit(f"ERROR: Missing required columns in variants CSV: {missing}")

    # Standardize names used downstream
    v = v.rename(columns={"contig": "contig", "position": "Position"})
    if "variantID" in v.columns and "PTC_ID" not in v.columns:
        v = v.rename(columns={"variantID": "PTC_ID"})
    if "PTC_ID" not in v.columns:
        v["PTC_ID"] = v.index.astype(str)

    v["Chromosome"] = v["contig"].astype(str)
    v["Position"]   = v["Position"].astype(int)
    v = harmonize_chr(v, "Chromosome", want_chr_prefix=True)

    # explode by txnames (strip version); handle multi-valued txnames
    rows = []
    for _, r in v.iterrows():
        txs = split_txnames(r["txnames"])
        if not txs:
            rows.append([r["PTC_ID"], r["Chromosome"], r["Position"], np.nan])
        else:
            for t in txs:
                rows.append([r["PTC_ID"], r["Chromosome"], r["Position"], strip_version(t)])
    v_exp = pd.DataFrame(rows, columns=["PTC_ID","Chromosome","Position","TxName"])
    return v[["PTC_ID","Chromosome","Position"]].drop_duplicates(), v_exp

def load_ejc(ejc_path):
    """
    EJC intervals as TSV/BED: chr  start  end  ...
    Only the first 3 columns are used.
    """
    ejc = pd.read_csv(ejc_path, sep="\t", header=None, comment="#")
    if ejc.shape[1] < 3:
        raise SystemExit("EJC file must have at least 3 columns: chr, start, end")
    ejc = ejc.iloc[:, :3]
    ejc.columns = ["Chromosome","Start","End"]
    ejc["Start"] = ejc["Start"].astype(int)
    ejc["End"]   = ejc["End"].astype(int)
    ejc = ejc[ejc["End"] > ejc["Start"]].copy()
    ejc = harmonize_chr(ejc, "Chromosome", want_chr_prefix=True)
    return pr.PyRanges(ejc)

def build_junctions_from_gtf(gtf_path):
    """
    Build exon–exon junctions (1-bp points) per transcript from GTF.
    transcript_id versions are stripped to match txnames provided.
    """
    gtf = pr.read_gtf(gtf_path)
    gdf = gtf.df.copy()
    if "transcript_id" not in gdf.columns:
        raise SystemExit("GTF is missing transcript_id attribute.")
    gdf["TranscriptID"] = gdf["transcript_id"].astype(str).map(strip_version)

    exons = gdf[gdf.Feature == "exon"][["Chromosome","Start","End","Strand","TranscriptID","exon_number"]].copy()

    # derive exon rank if needed
    if "exon_number" not in exons.columns or exons["exon_number"].isna().all():
        exons = (exons
                 .sort_values(["TranscriptID","Chromosome","Start","End"])
                 .assign(exon_rank=lambda d: d.groupby("TranscriptID").cumcount()+1))
    else:
        exons["exon_rank"] = (exons["exon_number"].astype(str)
                              .str.extract(r"(\d+)").fillna("1").astype(int))

    # create junctions as 1-bp [pos, pos+1)
    out = []
    for tx, d in exons.groupby("TranscriptID"):
        d = d.sort_values("exon_rank")
        if len(d) < 2:
            continue
        strand = d["Strand"].iloc[0]
        if strand == "+":
            jpos = d["End"].iloc[:-1].values
        else:
            jpos = d["Start"].iloc[:-1].values
        out.append(pd.DataFrame({
            "TranscriptID": tx,
            "Chromosome": d["Chromosome"].iloc[:-1].values,
            "Strand": strand,
            "Start": jpos,
            "End": jpos + 1
        }))
    junc = pd.concat(out, ignore_index=True) if out else pd.DataFrame(
        columns=["TranscriptID","Chromosome","Strand","Start","End"]
    )
    junc = harmonize_chr(junc, "Chromosome", want_chr_prefix=True)
    return junc

# ------------------ build windows around -24bp from nearest junction ------------------

def build_junction_minus24_windows(v_exp, junctions, offset=24, window_halfwidth=15):
    """
    For each exploded variant row (PTC_ID, Chromosome, Position, TxName),
    find the nearest downstream junction and build a window centered at -24bp from it:
      +: nearest_junction_pos - 24, then ±15bp window
      -: nearest_junction_pos + 24, then ±15bp window
    Returns PyRanges with: PTC_ID, Chromosome, Start, End, Strand, junction_pos
    """
    records = []
    j_by_tx = {k: v for k, v in junctions.groupby("TranscriptID")}
    tx_strand = (junctions[["TranscriptID","Strand"]]
                 .drop_duplicates()
                 .set_index("TranscriptID")["Strand"].to_dict())

    for _, r in v_exp.iterrows():
        tx = r["TxName"]
        if pd.isna(tx):
            continue  # cannot build tx-aware window
        chr_ = r["Chromosome"]
        pos  = int(r["Position"])
        pid  = r["PTC_ID"]

        js = j_by_tx.get(tx, None)
        strand = tx_strand.get(tx, None)
        if js is None or js.empty or strand is None:
            continue

        if strand == "+":
            # downstream junctions are those with Start > pos
            downstream = js[(js["Chromosome"] == chr_) & (js["Start"] > pos)]
            if downstream.empty:
                continue
            # nearest downstream junction
            nearest_junc = int(downstream["Start"].min())
            # -24bp from junction (upstream in genomic coords)
            target_pos = nearest_junc - offset
            # window: [target_pos - 15, target_pos + 15 + 1) for half-open interval
            start_win = target_pos - window_halfwidth
            end_win = target_pos + window_halfwidth + 1
            records.append([pid, chr_, start_win, end_win, strand, nearest_junc])
        else:  # negative strand
            # downstream junctions are those with Start < pos
            downstream = js[(js["Chromosome"] == chr_) & (js["Start"] < pos)]
            if downstream.empty:
                continue
            # nearest downstream junction
            nearest_junc = int(downstream["Start"].max())
            # -24bp from junction (upstream in transcription direction = downstream in genomic coords)
            target_pos = nearest_junc + offset
            # window: [target_pos - 15, target_pos + 15 + 1)
            start_win = target_pos - window_halfwidth
            end_win = target_pos + window_halfwidth + 1
            records.append([pid, chr_, start_win, end_win, strand, nearest_junc])

    cols = ["PTC_ID","Chromosome","Start","End","Strand","nearest_junction"]
    win_df = pd.DataFrame(records, columns=cols) if records else pd.DataFrame(columns=cols)
    return pr.PyRanges(win_df)

# ------------------ parallel overlap counting ------------------

def _count_on_chrom(args):
    chrom, win_df, ejc_df = args
    if win_df.empty or ejc_df.empty:
        return pd.DataFrame(columns=["PTC_ID","Chromosome","Start","End","Strand","nearest_junction","ejc_count_in_window","has_ejc_overlap"])
    gr_win = pr.PyRanges(win_df)
    gr_ejc = pr.PyRanges(ejc_df)
    counted = gr_win.count_overlaps(gr_ejc).df.rename(columns={"NumberOverlaps":"ejc_count_in_window"})
    # Add binary indicator: "yes" if count > 0, "no" otherwise
    counted["has_ejc_overlap"] = counted["ejc_count_in_window"].apply(lambda x: "yes" if x > 0 else "no")
    return counted[["PTC_ID","Chromosome","Start","End","Strand","nearest_junction","ejc_count_in_window","has_ejc_overlap"]]

def parallel_count_overlaps(wins: pr.PyRanges, ejc: pr.PyRanges, threads: int = 4):
    """Shard by chromosome and count overlaps in parallel workers."""
    if len(wins) == 0:
        return pd.DataFrame(columns=["PTC_ID","Chromosome","Start","End","Strand","nearest_junction","ejc_count_in_window","has_ejc_overlap"])
    wdf = wins.df
    edf = ejc.df
    chromosomes = sorted(set(wdf["Chromosome"]) & set(edf["Chromosome"]))
    tasks = []
    for ch in chromosomes:
        tasks.append((ch,
                      wdf[wdf["Chromosome"] == ch].copy(),
                      edf[edf["Chromosome"] == ch].copy()))
    out = []
    with ProcessPoolExecutor(max_workers=max(1, threads)) as ex:
        futures = [ex.submit(_count_on_chrom, t) for t in tasks]
        for f in as_completed(futures):
            out.append(f.result())
    if not out:
        return pd.DataFrame(columns=["PTC_ID","Chromosome","Start","End","Strand","nearest_junction","ejc_count_in_window","has_ejc_overlap"])
    return pd.concat(out, ignore_index=True)

# ------------------ main ------------------

def main():
    ap = argparse.ArgumentParser(
        description="EJC occupancy in 30bp window centered at -24bp from nearest downstream junction (CSV variants). Parallel by chromosome."
    )
    ap.add_argument("--variants", required=True, help="Variants CSV with columns: contig, position, txnames, [variantID]")
    ap.add_argument("--ejc", required=True, help="EJC intervals TSV/BED: chr start end ...")
    ap.add_argument("--gtf", required=True, help="Gencode v26 GTF (gz or not)")
    ap.add_argument("--out", required=True, help="Output TSV[.gz]")
    ap.add_argument("--offset", type=int, default=24, help="Distance upstream from junction (default 24 bp)")
    ap.add_argument("--window", type=int, default=15, help="Half-width of window around target position (default 15 bp)")
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 4, help="Parallel workers (default: all cores)")
    ap.add_argument("--aggregate", choices=["none","max","sum"], default="none",
                    help="Aggregate across multiple txnames per PTC (default: none)")
    args = ap.parse_args()

    print("[*] Loading variants CSV (columns: contig, position, txnames)…")
    v_orig, v_exp = load_variants_csv(args.variants)

    print("[*] Loading EJC intervals…")
    ejc = load_ejc(args.ejc)

    print("[*] Building exon–exon junctions from GTF…")
    junctions = build_junctions_from_gtf(args.gtf)

    print(f"[*] Building windows at -{args.offset}bp from nearest downstream junction (±{args.window}bp)…")
    wins = build_junction_minus24_windows(v_exp, junctions, offset=args.offset, window_halfwidth=args.window)

    ptc_with_win = set(wins.df["PTC_ID"].unique()) if len(wins) else set()
    all_ptc_ids  = set(v_orig["PTC_ID"])

    print(f"[*] Counting overlaps in parallel (threads={args.threads})…")
    occ_df = parallel_count_overlaps(wins, ejc, threads=args.threads)

    # Optional aggregation across txnames
    if args.aggregate != "none" and len(occ_df):
        if args.aggregate == "max":
            agg_df = occ_df.groupby("PTC_ID", as_index=False).agg({
                "ejc_count_in_window": "max"
            })
            # Set has_ejc_overlap based on aggregated count
            agg_df["has_ejc_overlap"] = agg_df["ejc_count_in_window"].apply(lambda x: "yes" if x > 0 else "no")
        else:  # sum
            agg_df = occ_df.groupby("PTC_ID", as_index=False).agg({
                "ejc_count_in_window": "sum"
            })
            # Set has_ejc_overlap based on aggregated count
            agg_df["has_ejc_overlap"] = agg_df["ejc_count_in_window"].apply(lambda x: "yes" if x > 0 else "no")
        occ_df = agg_df  # one row per PTC
    else:
        # keep per-(PTC, transcript window) rows with coordinates
        pass

    # Add zeros for PTCs that never got a window (likely last exon across all txnames)
    zero_ids = sorted(all_ptc_ids - set(occ_df["PTC_ID"].unique()))
    if zero_ids:
        zeros = v_orig.loc[v_orig["PTC_ID"].isin(zero_ids), ["PTC_ID","Chromosome"]].drop_duplicates().copy()
        if args.aggregate == "none":
            zeros["Start"] = np.nan
            zeros["End"] = np.nan
            zeros["Strand"] = np.nan
            zeros["nearest_junction"] = np.nan
        zeros["ejc_count_in_window"] = 0
        zeros["has_ejc_overlap"] = "no"
        occ_df = pd.concat([occ_df, zeros], ignore_index=True)

    # Sort & write
    if "Start" in occ_df.columns:
        occ_df = occ_df.sort_values(["PTC_ID","Chromosome","Start"], na_position="last")
    else:
        occ_df = occ_df.sort_values(["PTC_ID"])

    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    if outp.suffix.endswith("gz"):
        occ_df.to_csv(outp, sep="\t", index=False, compression="gzip")
    else:
        occ_df.to_csv(outp, sep="\t", index=False)

    print(f"[*] Done. Wrote {len(occ_df)} rows to {outp}")

if __name__ == "__main__":
    main()