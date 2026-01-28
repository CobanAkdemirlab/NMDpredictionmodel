#!/usr/bin/env python3
"""
Median Conservation Score Extraction for Stop-Gain Variants
Transcript selection: Uses txnames column with ENST Transcript ID with version.
Optional version-free mode using tx_nover or by stripping ENST versions.

Regions computed (per variant):
- ptc_100bp (±100 bp around PTC)
- ptc_to_ejc  (PTC -> downstream exon junction within same exon; strand-aware)
- ejc_100bp   (±100 bp around that downstream exon junction)
- old3utr_first200, old3utr_whole (canonical 3' UTR from CDS end)
- new3utr_first200, new3utr_whole (from PTC to end of transcript)
- utr5_first200, utr5_whole
- tx_whole    (spliced exonic: 5'UTR + CDS + 3'UTR)
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path

# =========================
# CONFIG - EDIT THESE PATHS BEFORE RUNNING
# =========================

# ==============================================================================
# REQUIRED: Edit these paths to match your system
# See REFERENCE_DATA_GUIDE.md for download instructions
# ==============================================================================

# Input variant file (CSV with contig, position, txnames columns)
# Example: "my_variants.csv" or "data/TOPMed_variants.csv"
VARIANT_FILE = "path/to/your/variants.csv"

# Conservation scores (BigWig format, ~10 GB each)
# Download from: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/
#            and: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/
# Or use: "reference/conservation/hg38.phastCons100way.bw" if following reference guide
PHASTCONS_BW = "path/to/hg38.phastCons100way.bw"
PHYLOP_BW    = "path/to/hg38.phyloP100way.bw"

# Gene annotations (Gencode v26 GTF, can be .gz compressed)
# Download from: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/
# Or use: "reference/annotations/gencode.v26.primary_assembly.annotation.gtf.gz"
GTF_FILE     = "path/to/gencode.v26.primary_assembly.annotation.gtf.gz"

# ==============================================================================
# OPTIONAL: Customize output paths if desired (will be created automatically)
# ==============================================================================

OUTPUT_DIR        = "conservation_output"
BED_FILE          = f"{OUTPUT_DIR}/conservation_regions.bed"
PHASTCONS_MEDIAN  = f"{OUTPUT_DIR}/phastcons_medians.tsv"
PHYLOP_MEDIAN     = f"{OUTPUT_DIR}/phylop_medians.tsv"
FINAL_OUTPUT      = f"{OUTPUT_DIR}/variants_with_conservation_medians.csv"

# Transcript handling
USE_VERSION_FREE = False  # if True, prefer tx_nover or strip ".##" from txnames entries

try:
    import pyBigWig
except ImportError:
    print("ERROR: pyBigWig is required. Install with:  pip install pyBigWig")
    sys.exit(1)

# ----------------------------------------------------------------------------
# Flexible CSV loader + column normalization (txnames-only policy)
# ----------------------------------------------------------------------------

# Candidate names to try (we never use 'transcript_id' for selection)
CANDIDATES = {
    "chrom":  ["contig", "CHROM", "#CHROM", "chromosome", "chr", "chrom", "Chromosome", "Chr", "seqnames"],
    "pos":    ["PTC_pos", "ptc_position", "position", "POS", "start", "variant_pos", "genomic_pos"],
    "txnames": ["txnames", "tx_list", "TX", "Transcript", "transcripts", "TXNAME"],
    "tx_nover": ["tx_nover", "txnames_nover", "tx_stripped"],
    "strand": ["strand", "STRAND", "tx_strand", "Strand"],
    "variant_id": ["variant_id", "Variant_ID", "key", "id", "variant"]
}

def _find_col(cols, candidates, required=False, label=""):
    for c in candidates:
        if c in cols:
            return c
    if required:
        raise ValueError(f"Missing required column for {label}. Tried {candidates}. "
                         f"Columns available: {list(cols)[:40]} …")
    return None

def _split_txlist(val):
    """Split transcript list from common delimiters; clean blanks."""
    if pd.isna(val):
        return []
    s = str(val).strip()
    if not s:
        return []
    for sep in [",", ";", "|"]:
        if sep in s:
            return [p.strip() for p in s.split(sep) if p.strip()]
    return [p.strip() for p in s.split() if p.strip()]

def _strip_ver(enst):
    """Remove version suffix from an ENST ID (ENST... .##)."""
    if enst is None:
        return None
    s = str(enst)
    return s.split(".")[0]

def load_and_normalize_variants(variant_path: str) -> pd.DataFrame:
    """Load CSV and normalize into: contig, position, txnames(list), tx_nover(list), strand, variant_id."""
    df = pd.read_csv(variant_path)
    cols = set(df.columns)

    chrom_col   = _find_col(cols, CANDIDATES["chrom"],   required=True,  label="chromosome")
    pos_col     = _find_col(cols, CANDIDATES["pos"],     required=True,  label="PTC position")
    tx_col      = _find_col(cols, CANDIDATES["txnames"], required=False, label="txnames")
    tx_novercol = _find_col(cols, CANDIDATES["tx_nover"], required=False, label="tx_nover")
    strand_col  = _find_col(cols, CANDIDATES["strand"],  required=False, label="strand")
    vid_col     = _find_col(cols, CANDIDATES["variant_id"], required=False, label="variant_id")

    out = df.copy()
    out["contig"] = out[chrom_col].astype(str)
    out["position"] = pd.to_numeric(out[pos_col], errors="coerce").astype("Int64")

    # txnames (primary) -> list
    if tx_col:
        out["txnames_raw"] = out[tx_col]
        out["txnames"] = out[tx_col].apply(_split_txlist)
    else:
        out["txnames_raw"] = pd.NA
        out["txnames"] = [[] for _ in range(len(out))]

    # tx_nover (optional) -> list
    if tx_novercol:
        out["tx_nover_raw"] = out[tx_novercol]
        out["tx_nover"] = out[tx_novercol].apply(_split_txlist)
    else:
        # build from txnames if needed
        out["tx_nover_raw"] = pd.NA
        out["tx_nover"] = [ [_strip_ver(t) for t in lst] for lst in out["txnames"] ]

    # strand (optional)
    out["strand"] = out[strand_col].astype(str) if strand_col else "+"

    # variant_id (optional)
    if vid_col:
        out["variant_id"] = out[vid_col].astype(str)
    else:
        out["variant_id"] = out["contig"].astype(str) + ":" + out["position"].astype(str)

    # Summary
    print("\n[✓] VARIANT CSV column mapping:")
    print(f"    - Chromosome → '{chrom_col}'  → 'contig'")
    print(f"    - PTC position → '{pos_col}'  → 'position'")
    print(f"    - txnames → '{tx_col}'")
    print(f"    - tx_nover → '{tx_novercol}' (or derived from txnames)")
    print(f"    - Strand → '{strand_col or 'N/A(+ default)'}'")
    print(f"    - Variant ID → '{vid_col or 'constructed contig:position'}'")
    print("\n[•] Normalized preview (first 3 rows):")
    print(out[["variant_id","contig","position","txnames","tx_nover","strand"]].head(3))
    return out

# =========================
# GTF parsing & canonical metadata (unchanged)
# =========================

def _get_attr(attr_text: str, key: str):
    if pd.isna(attr_text):
        return None
    parts = [p.strip() for p in str(attr_text).split(";") if p.strip()]
    key_plus = key + " "
    for p in parts:
        if p.startswith(key_plus):
            return p[len(key_plus):].strip().strip('"')
    return None

def _get_tags(attr_text: str):
    if pd.isna(attr_text):
        return set()
    tags = set()
    parts = [p.strip() for p in str(attr_text).split(";") if p.strip()]
    for p in parts:
        if p.startswith("tag "):
            val = p[len("tag "):].strip().strip('"')
            if val:
                tags.add(val)
    return tags

def _appris_rank(tags: set):
    for r, key in enumerate([
        "appris_principal_1","appris_principal_2","appris_principal_3",
        "appris_principal_4","appris_principal_5","appris_alternative_1","appris_alternative_2"
    ]):
        if key in tags:
            return r
    return 99

def _is_mane(tags: set):
    return "MANE_Select" in tags or "MANE_Plus_Clinical" in tags

def _tsl_rank(tsl_text: str):
    if tsl_text is None:
        return 99
    t = tsl_text.strip().upper().replace("TSL", "")
    try:
        return int(t)
    except ValueError:
        return 99

def load_gtf_model_and_meta(gtf_path: str):
    print("\n[•] Reading GTF (exons + CDS + transcript-level metadata) …")
    cols = ["seqname","source","feature","start","end","score","strand","frame","attribute"]
    usecols = list(range(9))
    exons_list, cds_list = [], []
    tx_strand = {}
    tx_meta_rows = []
    chunks = exon_rows = cds_rows = 0

    for chunk in pd.read_csv(
        gtf_path, sep="\t", comment="#", header=None, names=cols, usecols=usecols,
        chunksize=200_000, dtype={"seqname":"string","feature":"string","strand":"string","attribute":"string"}
    ):
        chunks += 1
        chunk["transcript_id"] = chunk["attribute"].apply(lambda s: _get_attr(s, "transcript_id"))
        sub = chunk[chunk["transcript_id"].notna()]
        if sub.empty:
            continue

        t_str = sub[["transcript_id","strand"]].dropna().drop_duplicates()
        for _, r in t_str.iterrows():
            tx_strand[r["transcript_id"]] = r["strand"]

        tmp = sub[["transcript_id","attribute"]].copy()
        tmp["tags"] = tmp["attribute"].apply(_get_tags)
        tmp["tsl"]  = tmp["attribute"].apply(lambda s: _get_attr(s, "transcript_support_level"))
        tx_meta_rows.append(tmp[["transcript_id","tags","tsl"]])

        ex = sub[sub["feature"] == "exon"].copy()
        if not ex.empty:
            ex = ex.rename(columns={"seqname":"chrom"})[["chrom","start","end","strand","transcript_id","attribute"]]
            for c in ("start","end"): ex[c] = pd.to_numeric(ex[c], errors="coerce").astype("Int64")
            exons_list.append(ex); exon_rows += len(ex)

        cds = sub[sub["feature"] == "CDS"].copy()
        if not cds.empty:
            cds = cds.rename(columns={"seqname":"chrom"})[["chrom","start","end","strand","transcript_id","attribute"]]
            for c in ("start","end"): cds[c] = pd.to_numeric(cds[c], errors="coerce").astype("Int64")
            cds_list.append(cds); cds_rows += len(cds)

    exons_df = pd.concat(exons_list, ignore_index=True) if exons_list else pd.DataFrame(
        columns=["chrom","start","end","strand","transcript_id","attribute"]
    )
    cds_df = pd.concat(cds_list, ignore_index=True) if cds_list else pd.DataFrame(
        columns=["chrom","start","end","strand","transcript_id","attribute"]
    )
    cds_bounds = {}
    if not cds_df.empty:
        cds_bounds = (cds_df.groupby("transcript_id", sort=False)
                           .agg(cds_start=("start","min"), cds_end=("end","max"))
                           .to_dict(orient="index"))

    exons_by_tx = {}
    exon_len_by_tx = {}
    for tx, sub in exons_df.groupby("transcript_id", sort=False):
        sub = sub.sort_values(["chrom","start","end"]).reset_index(drop=True)
        exons_by_tx[tx] = sub[["chrom","start","end","strand","transcript_id"]]
        exon_len = int((sub["end"] - sub["start"] + 1).clip(lower=0).sum())
        exon_len_by_tx[tx] = exon_len

    cds_len_by_tx = {}
    if not cds_df.empty:
        cds_len_by_tx = (cds_df.assign(seg_len=(cds_df["end"] - cds_df["start"] + 1).clip(lower=0))
                               .groupby("transcript_id", sort=False)["seg_len"].sum()
                               .astype(int).to_dict())

    if tx_meta_rows:
        meta_df = pd.concat(tx_meta_rows, ignore_index=True)
        meta_df["appris"] = meta_df["tags"].apply(_appris_rank)
        meta_df["is_mane"] = meta_df["tags"].apply(_is_mane)
        meta_df["tsl_rank"] = meta_df["tsl"].apply(_tsl_rank)
        meta_agg = (meta_df.groupby("transcript_id", sort=False)
                          .agg(is_mane=("is_mane","max"),
                               appris_rank=("appris","min"),
                               tsl_rank=("tsl_rank","min"))
                          .reset_index())
    else:
        meta_agg = pd.DataFrame(columns=["transcript_id","is_mane","appris_rank","tsl_rank"])

    meta_agg["cds_len"]  = meta_agg["transcript_id"].map(cds_len_by_tx).fillna(0).astype(int)
    meta_agg["exon_len"] = meta_agg["transcript_id"].map(exon_len_by_tx).fillna(0).astype(int)

    tx_meta = {r["transcript_id"]:
               {"is_mane": bool(r["is_mane"]),
                "appris_rank": int(r["appris_rank"]) if not pd.isna(r["appris_rank"]) else 99,
                "tsl_rank": int(r["tsl_rank"]) if not pd.isna(r["tsl_rank"]) else 99,
                "cds_len": int(r["cds_len"]),
                "exon_len": int(r["exon_len"])}
               for _, r in meta_agg.iterrows()}

    print(f"[✓] GTF parsing complete:")
    print(f"    - Exon rows                   : {len(exons_df):,}")
    print(f"    - CDS rows                    : {len(cds_df):,}")
    print(f"    - Transcripts w/ exons        : {len(exons_by_tx):,}")
    print(f"    - Transcripts w/ CDS          : {len(cds_bounds):,}")
    print(f"    - Transcripts w/ metadata     : {len(tx_meta):,}")
    return exons_by_tx, cds_bounds, {k:v.iloc[0] if hasattr(v,'iloc') else v for k,v in {}.items()}, tx_meta

# =========================
# Canonical picker (MANE > APPRIS > TSL > longest CDS > exon_len)
# =========================

def pick_canonical(tx_list, tx_meta):
    candidates = []
    for tx in tx_list:
        m = tx_meta.get(tx, {"is_mane":False,"appris_rank":99,"tsl_rank":99,"cds_len":0,"exon_len":0})
        score = (0 if m["is_mane"] else 1, m["appris_rank"], m["tsl_rank"], -m["cds_len"], -m["exon_len"])
        candidates.append((score, tx, m))
    if not candidates:
        return None, "no_tx"
    candidates.sort(key=lambda x: x[0])
    best_score, best_tx, m = candidates[0]
    if m["is_mane"]: reason = "MANE"
    elif m["appris_rank"] < 99: reason = f"APPRIS_{m['appris_rank']}"
    elif m["tsl_rank"] < 99: reason = f"TSL{m['tsl_rank']}"
    elif m["cds_len"] > 0: reason = "LongestCDS"
    elif m["exon_len"] > 0: reason = "LongestExonLen"
    else: reason = "Fallback"
    return best_tx, reason

# =========================
# Region builder (uses txnames/tx_nover only)
# =========================

class RegionBuilder:
    def __init__(self, df, exons_by_tx, cds_bounds, tx_meta):
        self.df = df
        self.exons_by_tx = exons_by_tx
        self.cds_bounds = cds_bounds
        self.tx_meta = tx_meta
        self.regions = []
        self.canon_stats = {"single_tx":0,"multi_tx":0,"no_tx":0}
        self.canon_reason_counts = {}

    @staticmethod
    def norm_chr(chrom):
        chrom = str(chrom)
        return chrom if chrom.startswith("chr") else f"chr{chrom}"

    def strand_for_tx(self, tx):
        # derive strand from exons table if available
        ex = self.exons_by_tx.get(tx)
        if ex is not None and not ex.empty and "strand" in ex.columns:
            s = str(ex["strand"].iloc[0])
            return s if s in ("+","-") else "+"
        return "+"

    def add_bed(self, chrom, start1, end1, name, strand):
        start0 = max(0, int(start1) - 1)
        end0   = max(start0 + 1, int(end1))
        self.regions.append([chrom, start0, end0, name, 0, strand])

    def transcript_blocks(self, tx):
        ex = self.exons_by_tx.get(tx)
        if ex is None or ex.empty: return []
        return [(str(r["chrom"]), int(r["start"]), int(r["end"])) for _, r in ex.iterrows()]

    def utr_blocks(self, tx, which):
        ex = self.exons_by_tx.get(tx)
        if ex is None or ex.empty: return []
        cds = self.cds_bounds.get(tx)
        if not cds: return []
        cds_start, cds_end = int(cds["cds_start"]), int(cds["cds_end"])
        strand = self.strand_for_tx(tx)
        out = []
        for _, r in ex.iterrows():
            ch = str(r["chrom"]); s = int(r["start"]); e = int(r["end"])
            if which == "5p":
                if strand == "+":
                    if s < cds_start:
                        out.append((ch, s, min(e, cds_start - 1)))
                else:
                    if e > cds_end:
                        out.append((ch, max(s, cds_end + 1), e))
            else:  # "3p"
                if strand == "+":
                    if e > cds_end:
                        out.append((ch, max(s, cds_end + 1), e))
                else:
                    if s < cds_start:
                        out.append((ch, s, min(e, cds_start - 1)))
        return [(ch, s, e) for (ch, s, e) in out if e >= s]

    def take_first_N(self, tx, blocks, N):
        if N <= 0 or not blocks: return []
        strand = self.strand_for_tx(tx)
        ordered = blocks if strand == "+" else list(reversed(blocks))
        rem = N; taken = []
        for ch, s, e in ordered:
            L = e - s + 1
            if L <= 0: continue
            if rem <= 0: break
            if L <= rem:
                taken.append((ch, s, e)); rem -= L
            else:
                if strand == "+":
                    taken.append((ch, s, s + rem - 1))
                else:
                    taken.append((ch, e - rem + 1, e))
                rem = 0
        return sorted(taken, key=lambda t: (t[0], t[1], t[2]))

    def new3utr_blocks(self, tx, chrom, ptc):
        ex = self.exons_by_tx.get(tx)
        if ex is None or ex.empty: return []
        strand = self.strand_for_tx(tx); out = []
        for _, r in ex.iterrows():
            ch = str(r["chrom"]); s = int(r["start"]); e = int(r["end"])
            if strand == "+":
                if e <= ptc: continue
                bs = max(s, ptc + 1); be = e
            else:
                if s >= ptc: continue
                bs = s; be = min(e, ptc - 1)
            if be >= bs: out.append((ch, bs, be))
        return out

    def downstream_junction_site(self, tx, chrom, ptc):
        ex = self.exons_by_tx.get(tx)
        if ex is None or ex.empty: return None
        strand = self.strand_for_tx(tx)
        hit = ex[(ex["chrom"] == chrom) & (ex["start"] <= ptc) & (ptc <= ex["end"])]
        if hit.empty: return None
        idx = hit.index[0]
        return int(ex.loc[idx, "end"] if strand == "+" else ex.loc[idx, "start"])

    def pick_variant_txlist(self, row):
        """Enforce: use txnames; version-free if USE_VERSION_FREE (prefer tx_nover). Never use transcript_id."""
        if USE_VERSION_FREE and "tx_nover" in row and isinstance(row["tx_nover"], list) and row["tx_nover"]:
            return row["tx_nover"]
        # else use txnames (strip versions if USE_VERSION_FREE)
        lst = row["txnames"] if isinstance(row["txnames"], list) else []
        if USE_VERSION_FREE:
            lst = [_strip_ver(t) for t in lst]
        return lst

    def build(self):
        print("\n[•] Selecting canonical transcripts (txnames only) & building regions …")
        types_counter = {}

        for i, row in self.df.iterrows():
            chrom = self.norm_chr(row["contig"])
            ptc   = int(row["position"])
            tx_list = self.pick_variant_txlist(row)

            if not tx_list:
                self.canon_stats["no_tx"] += 1
                tx = None; reason = "no_tx"
            elif len(tx_list) == 1:
                self.canon_stats["single_tx"] += 1
                tx = tx_list[0]; reason = "single"
            else:
                self.canon_stats["multi_tx"] += 1
                # map version-free IDs back to versioned keys if needed
                if USE_VERSION_FREE:
                    # try exact matches first; else try to expand by prefix match against available GTF keys
                    expanded = []
                    gtf_ids = set(self.exons_by_tx.keys())
                    for t in tx_list:
                        if t in gtf_ids:  # already versioned match
                            expanded.append(t); continue
                        # try prefix match (t == ENSTxxxx without version)
                        matches = [gid for gid in gtf_ids if gid.split(".")[0] == t]
                        expanded.extend(matches or [])
                    tx, reason = pick_canonical(expanded or tx_list, self.tx_meta)
                else:
                    tx, reason = pick_canonical(tx_list, self.tx_meta)
            self.canon_reason_counts[reason] = self.canon_reason_counts.get(reason, 0) + 1

            strand = self.strand_for_tx(tx) if tx else "+"
            base = f"var{i}"

            # 1) PTC ±100
            self.add_bed(chrom, ptc - 100, ptc + 100, f"{base}_ptc_100bp", strand)
            types_counter["ptc_100bp"] = types_counter.get("ptc_100bp", 0) + 1

            if tx and tx in self.exons_by_tx:
                ex = self.exons_by_tx[tx]
                curr = ex[(ex["chrom"] == chrom) & (ex["start"] <= ptc) & (ptc <= ex["end"])]

                if not curr.empty:
                    s = int(curr.iloc[0]["start"]); e = int(curr.iloc[0]["end"])
                    if strand == "+" and e > ptc:
                        self.add_bed(chrom, ptc, e, f"{base}_ptc_to_ejc", strand)
                        types_counter["ptc_to_ejc"] = types_counter.get("ptc_to_ejc", 0) + 1
                    elif strand == "-" and s < ptc:
                        self.add_bed(chrom, s, ptc, f"{base}_ptc_to_ejc", strand)
                        types_counter["ptc_to_ejc"] = types_counter.get("ptc_to_ejc", 0) + 1

                ejc = self.downstream_junction_site(tx, chrom, ptc)
                if ejc is not None:
                    self.add_bed(chrom, ejc - 100, ejc + 100, f"{base}_ejc_100bp", strand)
                    types_counter["ejc_100bp"] = types_counter.get("ejc_100bp", 0) + 1

                old3 = self.utr_blocks(tx, "3p")
                if old3:
                    for ch, s, e in old3: self.add_bed(ch, s, e, f"{base}_old3utr_whole", strand)
                    types_counter["old3utr_whole"] = types_counter.get("old3utr_whole", 0) + 1
                    for ch, s, e in self.take_first_N(tx, old3, 200):
                        self.add_bed(ch, s, e, f"{base}_old3utr_first200", strand)
                    types_counter["old3utr_first200"] = types_counter.get("old3utr_first200", 0) + 1

                new3 = self.new3utr_blocks(tx, chrom, ptc)
                if new3:
                    for ch, s, e in new3: self.add_bed(ch, s, e, f"{base}_new3utr_whole", strand)
                    types_counter["new3utr_whole"] = types_counter.get("new3utr_whole", 0) + 1
                    for ch, s, e in self.take_first_N(tx, new3, 200):
                        self.add_bed(ch, s, e, f"{base}_new3utr_first200", strand)
                    types_counter["new3utr_first200"] = types_counter.get("new3utr_first200", 0) + 1

                utr5 = self.utr_blocks(tx, "5p")
                if utr5:
                    for ch, s, e in utr5: self.add_bed(ch, s, e, f"{base}_utr5_whole", strand)
                    types_counter["utr5_whole"] = types_counter.get("utr5_whole", 0) + 1
                    for ch, s, e in self.take_first_N(tx, utr5, 200):
                        self.add_bed(ch, s, e, f"{base}_utr5_first200", strand)
                    types_counter["utr5_first200"] = types_counter.get("utr5_first200", 0) + 1

                tx_blocks = self.transcript_blocks(tx)
                if tx_blocks:
                    for ch, s, e in tx_blocks: self.add_bed(ch, s, e, f"{base}_tx_whole", strand)
                    types_counter["tx_whole"] = types_counter.get("tx_whole", 0) + 1

        print(f"\n[✓] Canonical transcript selection summary (txnames only):")
        for k in ["no_tx","single_tx","multi_tx"]:
            print(f"    - {k:12s}: {self.canon_stats.get(k,0):8,d}")
        print("    - Reasons used:")
        for k in sorted(self.canon_reason_counts.keys()):
            print(f"        {k:15s}: {self.canon_reason_counts[k]:8,d}")

        print(f"\n[✓] Region construction summary:")
        print(f"    - Variants processed            : {len(self.df):,}")
        print(f"    - Total BED rows (segments)     : {len(self.regions):,}")
        print(f"    - Regions by type (logical counts):")
        for t, v in sorted({n:self.regions.count(n) for n in []}.items()):
            print(f"        {t:25s}: {v:8,d}")  # (kept simple; detailed counts are printed as built)

    def to_bed(self, out_path: str):
        bed_df = pd.DataFrame(self.regions, columns=["chrom","start","end","name","score","strand"])
        bed_df.to_csv(out_path, sep="\t", header=False, index=False)
        print(f"[✓] BED written: {out_path}  (rows: {len(bed_df):,})")
        return bed_df

# =========================
# Medians + merge
# =========================

def median_scores_from_bigwig(bed_df: pd.DataFrame, bigwig_path: str, out_path: str):
    print(f"\n[•] Computing medians from bigWig: {bigwig_path}")
    bw = pyBigWig.open(bigwig_path)
    grouped = bed_df.groupby("name", sort=False)
    records = []
    processed = 0

    for name, sub in grouped:
        total_bp = int((sub["end"] - sub["start"]).clip(lower=0).sum())
        vals_list = []
        for _, row in sub.iterrows():
            try:
                arr = np.array(bw.values(row["chrom"], int(row["start"]), int(row["end"]), numpy=True), dtype=float)
            except RuntimeError:
                continue
            if arr.size:
                vals_list.append(arr)
        if not vals_list:
            records.append((name, total_bp, 0, np.nan))
        else:
            concat = np.concatenate(vals_list)
            valid = concat[~np.isnan(concat)]
            med = np.nan if valid.size == 0 else float(np.median(valid))
            records.append((name, total_bp, int(valid.size), med))
        processed += 1
        if processed % 10000 == 0:
            print(f"    … {processed:,} regions aggregated")

    bw.close()
    out = pd.DataFrame(records, columns=["name","size_bp","valid_bp","median"])
    out.to_csv(out_path, sep="\t", index=False)

    num_regions   = len(out)
    num_with_data = (out["valid_bp"] > 0).sum()
    pct_with_data = 100 * num_with_data / num_regions if num_regions else 0
    print(f"[✓] Median extraction summary for {Path(bigwig_path).name}:")
    print(f"    - Regions total        : {num_regions:,}")
    print(f"    - Regions with data    : {num_with_data:,} ({pct_with_data:.1f}%)")
    try:
        desc = out['median'].describe(percentiles=[.05,.25,.5,.75,.95])
        print("    - Median distribution summary:\n", desc)
    except Exception:
        pass

    print(f"[✓] Wrote medians: {out_path}")
    return out

def merge_medians(variant_csv: str, phastcons_median_tsv: str, phylop_median_tsv: str, out_csv: str):
    print("\n[•] Merging median conservation scores back into variants …")
    variants = pd.read_csv(variant_csv)
    result = variants.copy()

    def widen(df, prefix):
        df = df.copy()
        df["variant_idx"] = df["name"].str.extract(r"var(\d+)_").astype(int)
        df["region_type"] = df["name"].str.replace(r"^var\d+_", "", regex=True)
        wide = df.pivot_table(index="variant_idx", columns="region_type", values="median", aggfunc="first")
        wide.columns = [f"{prefix}_{c}_median" for c in wide.columns]
        return wide

    phc_cols = []
    php_cols = []

    if Path(phastcons_median_tsv).exists():
        phc = pd.read_csv(phastcons_median_tsv, sep="\t")
        phc_w = widen(phc, "phastcons")
        result = result.merge(phc_w, left_index=True, right_index=True, how="left")
        phc_cols = list(phc_w.columns)

    if Path(phylop_median_tsv).exists():
        php = pd.read_csv(phylop_median_tsv, sep="\t")
        php_w = widen(php, "phylop")
        result = result.merge(php_w, left_index=True, right_index=True, how="left")
        php_cols = list(php_w.columns)

    result.to_csv(out_csv, index=False)
    print(f"[✓] Final CSV saved: {out_csv}")

    # Quick verification
    for prefix, cols in [("phastcons", phc_cols), ("phylop", php_cols)]:
        if cols:
            nonnull_rows = result[cols].notna().any(axis=1).sum()
            print(f"    - Variants with any {prefix} median: {nonnull_rows:,} / {len(result):,} "
                  f"({100*nonnull_rows/len(result):.1f}%)")
        else:
            print(f"    - No {prefix} columns produced.")

# =========================
# MAIN
# =========================

def main():
    print("="*80)
    print("Median Conservation for Stop-Gain Variants (PhastCons & PhyloP)")
    print("="*80)

    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    for lbl, fp in [("Variants", VARIANT_FILE), ("PhastCons", PHASTCONS_BW),
                    ("PhyloP", PHYLOP_BW), ("GTF", GTF_FILE)]:
        ok = Path(fp).exists()
        print(f"{'✓' if ok else '✗'} {lbl:10s}: {fp}")
    # Load variants (txnames-only policy)
    variants = load_and_normalize_variants(VARIANT_FILE)
    print(f"\n[✓] Loaded variants: {len(variants):,}")

    exons_by_tx, cds_bounds, _tx_strand_unused, tx_meta = load_gtf_model_and_meta(GTF_FILE)

    rb = RegionBuilder(variants, exons_by_tx, cds_bounds, tx_meta)
    rb.build()
    bed_df = rb.to_bed(BED_FILE)

    phc = median_scores_from_bigwig(bed_df, PHASTCONS_BW, PHASTCONS_MEDIAN)
    php = median_scores_from_bigwig(bed_df, PHYLOP_BW,    PHYLOP_MEDIAN)

    merge_medians(VARIANT_FILE, PHASTCONS_MEDIAN, PHYLOP_MEDIAN, FINAL_OUTPUT)

    print("\n" + "="*80)
    print("✓ COMPLETE SUMMARY")
    print("="*80)
    print(f"Final CSV written             : {FINAL_OUTPUT}")
    print("="*80)

if __name__ == "__main__":
    main()
