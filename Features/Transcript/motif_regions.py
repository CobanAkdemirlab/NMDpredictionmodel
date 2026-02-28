#!/usr/bin/env python3
"""

This script execute the following:
- Builds per-variant GRanges-like lists using gffutils:
  - CDS from PTC to end + all 3'UTR
  - True 3'UTR extracted from the annotation (or exon–CDS arithmetic fallback)
  - First 200 bp of 3'UTR (strand-aware)
  - "newUTR": CDS segment from PTC to first downstream exon–exon junction
  - First 200 bp of newUTR (strand-aware)
  - ±100 bp windows around PTC and the first downstream EJC
  - PTC → EJC segment
- Writes FIMO batches with >seq_{row_index} headers so sequence_name maps exactly back to rows.
- Maps motifs to genes using ATtRACT Matrix_id; falls back to parsing motif_id like "<id>:<GENE>:..."
- DOES NOT modify the variants CSV or its row order/content in any way; outputs are separate.

Outputs:
  - fimo_regions/<region>/batch_*/results/fimo.tsv
  - matrices/motif_presence.csv (0/1 per <region>.<Gene> per variant row)
"""

import os, re, time, shutil, subprocess
from pathlib import Path
from multiprocessing import Pool, cpu_count

import pandas as pd
import numpy as np
import gffutils
from Bio.Seq import Seq
from pyfaidx import Fasta

# ---------------- Config (EDIT paths as needed) ----------------
WORKDIR   = Path(".").resolve()
VARIANTS  = WORKDIR / "/Users/iegab/TOPMed_stopgain_final_Feb132026.csv"
GTF       = WORKDIR / "/Users/iegab/Downloads/gencode.v26.primary_assembly.annotation.gtf.gz"
GENOME    = WORKDIR / "/Users/iegab/Downloads/hg38.fa"
MOTIFS    = WORKDIR / "/Users/iegab/Downloads/memefile.txt"                 # MEME format
ATTRACT   = WORKDIR / "/Users/iegab/Downloads/ATtRACT_db_human.txt"         # tab-delimited
MEME_BIN  = None                                      # e.g., "/usr/local/meme/bin" or None
IS_RNA    = False                                     # True if your motif file expects U instead of T

FIMO_THRESHOLD   = "1e-4"
FIMO_TIMEOUT_S   = 1200
FIMO_BATCH_SIZE  = 400
N_CORES          = max(1, cpu_count() - 1)

# ---- MEME PATH prepend like the user's system layout ----
try:
    from pathlib import Path as _P
    _meme_bin = str(_P.home() / "meme" / "bin")
    _meme_libexec = str(_P.home() / "meme" / "libexec" / "meme-5.5.5")
    for _d in (_meme_bin, _meme_libexec):
        if os.path.isdir(_d) and _d not in os.environ.get("PATH",""):
            os.environ["PATH"] = f"{_d}:{os.environ['PATH']}"
except Exception:
    pass
if MEME_BIN and MEME_BIN not in os.environ.get("PATH",""):
    os.environ["PATH"] = f"{MEME_BIN}:{os.environ['PATH']}"

# ---------------- Helpers ----------------
def _check_fimo():
    out = subprocess.run(["fimo", "--version"], capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError("FIMO not found on PATH. Adjust MEME_BIN or install MEME under ~/meme.")
    return (out.stdout.strip() or out.stderr.strip())

def _add_chr(chrom: str) -> str:
    return chrom if str(chrom).startswith("chr") else f"chr{chrom}"

def _get_tx(db, txname: str):
    # txnames may be "ENST...;ENST..." → take first
    txid = str(txname).split(";")[0].strip()
    return db[txid], txid

def _cds_list(db, tx):
    return list(db.children(tx, featuretype="CDS", order_by="start"))

def _exon_list(db, tx):
    return list(db.children(tx, featuretype="exon", order_by="start"))

def _utr3_from_schema(db, tx):
    """Try to fetch 3'UTR features using multiple common schemas."""
    for ft in ("three_prime_UTR", "three_prime_utr", "3UTR"):
        feats = list(db.children(tx, featuretype=ft, order_by="start"))
        if feats:
            return feats
    # Some GTFs use generic UTR with attributes signaling 3'UTR.
    utrs = list(db.children(tx, featuretype="UTR", order_by="start"))
    three_like = []
    for u in utrs:
        attrs = getattr(u, "attributes", {}) or {}
        labels = []
        for k in ("utr_type", "UTR_type", "type", "tag", "biotype"):
            v = attrs.get(k, [])
            if isinstance(v, (list, tuple)):
                labels.extend([str(x).lower() for x in v])
            elif v:
                labels.append(str(v).lower())
        label_str = " ".join(labels)
        if "three" in label_str and "utr" in label_str:
            three_like.append(u)
    return three_like

def _utr3_list(db, tx):
    feats = _utr3_from_schema(db, tx)
    if feats:
        return feats

    exons = _exon_list(db, tx)
    cds   = _cds_list(db, tx)
    if not exons or not cds:
        return []
    strand = tx.strand

    if strand == '+':
        last_cds_end = max(c.end for c in cds)
        out = []
        for e in exons:
            if e.end <= last_cds_end:
                continue
            s = max(e.start, last_cds_end + 1)
            if s <= e.end:
                part = type(e)(
                    seqid=e.seqid, source=e.source, featuretype="three_prime_UTR",
                    start=s, end=e.end, score=e.score, strand=e.strand,
                    frame=e.frame, attributes=e.attributes
                )
                out.append(part)
        return out
    else:
        first_cds_start = min(c.start for c in cds)
        out = []
        for e in exons:
            if e.start >= first_cds_start:
                continue
            e_end = min(e.end, first_cds_start - 1)
            if e.start <= e_end:
                part = type(e)(
                    seqid=e.seqid, source=e.source, featuretype="three_prime_UTR",
                    start=e.start, end=e_end, score=e.score, strand=e.strand,
                    frame=e.frame, attributes=e.attributes
                )
                out.append(part)
        return out

def _extract_seq(genome: Fasta, chrom: str, start: int, end: int, strand: str) -> str:
    if start > end:
        start, end = end, start
    seq = genome[chrom][start-1:end].seq
    if strand == '-':
        return str(Seq(seq).reverse_complement())
    return str(seq)

def _pieces_to_seq(genome: Fasta, pieces):
    if not pieces: return ""
    out = []
    for chrom, s, e, strand in pieces:
        try:
            out.append(_extract_seq(genome, chrom, s, e, strand))
        except Exception:
            return ""
    return "".join(out)

def _slice_firstN(gr_list, toward="start", bp=200):
    """Take first N bp across a list of regions following transcript direction."""
    if not gr_list: return []
    remaining = bp
    out = []
    it = gr_list if toward == "start" else list(reversed(gr_list))
    for region in it:
        seg_len = region.end - region.start + 1
        if remaining <= 0: break
        if seg_len <= remaining:
            out.append(region); remaining -= seg_len
        else:
            if toward == "start":
                part_start, part_end = region.start, region.start + remaining - 1
            else:
                part_start, part_end = region.end - remaining + 1, region.end
            part = type(region)(
                seqid=region.seqid, source=region.source, featuretype=region.featuretype,
                start=part_start, end=part_end, score=region.score, strand=region.strand,
                frame=region.frame, attributes=region.attributes
            )
            out.append(part); remaining = 0
    return out if toward == "start" else list(reversed(out))

def _first_ejc_downstream(db, tx, ptc_pos):
    """
    R-faithful EJC proxy:
    - '+' strand: EJC anchor = end of CDS block that contains the PTC
    - '-' strand: EJC anchor = start of CDS block that contains the PTC
    This matches your R: exon.cor[mut.exon] - coding.pos
    """
    cds = _cds_list(db, tx)
    if not cds:
        return None

    strand = tx.strand
    chrom = tx.chrom

    for c in cds:
        if c.start <= ptc_pos <= c.end:
            if strand == '+':
                return (chrom, c.end, strand)
            else:
                return (chrom, c.start, strand)

    return None

# ---------- REGION CONSTRUCTION ----------
def _build_regions_for_variant(genome, db, idx, chrom_raw, pos, txnames):
    out = {}
    chrom = _add_chr(str(chrom_raw))
    try:
        tx, txid = _get_tx(db, txnames)
    except Exception:
        return out
    strand = tx.strand
    # ---------- helpers local ----------
    def _cds_from_to_end(start_pos):
        """CDS pieces from start_pos (EJC) to CDS terminus, transcript-aware."""
        cds = _cds_list(db, tx)
        if not cds:
            return []
        pieces = []
        if strand == '+':
            for c in cds:
                if c.end < start_pos:
                    continue
                s = max(c.start, start_pos)
                e = c.end
                if s <= e:
                    pieces.append((chrom, s, e, strand))
        else:
            for c in cds:
                if c.start > start_pos:
                    continue
                s = c.start
                e = min(c.end, start_pos)
                if s <= e:
                    pieces.append((chrom, s, e, strand))
        return pieces

    def _mk_reg_like(template, s, e):
        return type(template)(
            seqid=tx.chrom, source="custom", featuretype="tmp",
            start=int(s), end=int(e), score=".", strand=strand,
            frame=".", attributes={}
        )
    # =========================================================
    # 1) 3'UTR (all + first200)
    # =========================================================
    utr3 = _utr3_list(db, tx)
    utr3_pieces = [(chrom, u.start, u.end, strand) for u in utr3] if utr3 else []
    if utr3_pieces:
        seq_utr3_all = _pieces_to_seq(genome, utr3_pieces)
        if seq_utr3_all:
            out["utr3_all"] = seq_utr3_all

        first200 = _slice_firstN(utr3, toward=("start" if strand == '+' else "end"), bp=200)
        if first200:
            seq_utr3_200 = _pieces_to_seq(genome, [(chrom, r.start, r.end, strand) for r in first200])
            if seq_utr3_200:
                out["utr3_200"] = seq_utr3_200

    # =========================================================
    # find EJC (required for ptc_to_ejc and new3utr)
    # =========================================================
    cds_list = _cds_list(db, tx)
    ejc = _first_ejc_downstream(db, tx, pos) if cds_list else None
    # DEBUG (only for first variant to avoid flooding output)
    if idx == 0 and ejc:
        _, jpos_debug, _ = ejc
        print("DEBUG VARIANT 0")
        print("PTC:", pos)
        print("EJC position returned:", jpos_debug)

    # =========================================================
    # 2) PTC → EJC 
    # =========================================================

    ptc2ejc_pieces = []

    if cds_list and ejc:
        _, jpos, _ = ejc
    
        if idx == 0:
            print("EJC used in ptc_to_ejc:", jpos)

        pieces = []

        for cds in cds_list:
            if strand == '+':
                if cds.end < pos:
                    continue

                s = max(cds.start, pos) if (cds.start <= pos <= cds.end) else cds.start
                e = cds.end

                if s <= jpos <= e:
                    e = jpos
                    pieces.append((chrom, s, e, strand))
                    break

                pieces.append((chrom, s, e, strand))

            else:
                if cds.start > pos:
                    continue

                s = cds.start
                e = min(cds.end, pos) if (cds.start <= pos <= cds.end) else cds.end

                if s <= jpos <= e:
                    s = jpos
                    pieces.append((chrom, s, e, strand))
                    break

                pieces.append((chrom, s, e, strand))

        ptc2ejc_pieces = pieces
        seq_ptc_to_ejc = _pieces_to_seq(genome, ptc2ejc_pieces)

        if seq_ptc_to_ejc:
            out["ptc_to_ejc"] = seq_ptc_to_ejc
    # =========================================================
    # 3) NEW 3'UTR (your corrected definition)
    #    = (PTC→EJC pieces) + (EJC→end of CDS) + (3'UTR)
    #    WITH 200bp version
    # =========================================================
    if ptc2ejc_pieces and ejc:
        _, jpos, _ = ejc
        ejc2cdsend_pieces = _cds_from_to_end(jpos)

        new3utr_pieces = ptc2ejc_pieces + ejc2cdsend_pieces + utr3_pieces
        seq_new3utr_all = _pieces_to_seq(genome, new3utr_pieces)
        if seq_new3utr_all:
            out["newutr_all"] = seq_new3utr_all

            # first 200 bp of new3utr in transcript direction
            template = cds_list[0] if cds_list else (utr3[0] if utr3 else None)
            if template is not None:
                reg_like = [_mk_reg_like(template, s, e) for (_, s, e, _) in new3utr_pieces]
                first200 = _slice_firstN(reg_like, toward=("start" if strand == '+' else "end"), bp=200)
                if first200:
                    seq_new3utr_200 = _pieces_to_seq(genome, [(chrom, r.start, r.end, strand) for r in first200])
                    if seq_new3utr_200:
                        out["newutr_200"] = seq_new3utr_200
    # =========================================================
    # 4a) PTC ±100 windows (genomic)
    # =========================================================
    
    out["ptc_pm100"] = _extract_seq(genome, chrom, pos-100, pos+99, strand)

    # =========================================================
    # 4b) EJC ±100 window (genomic, symmetric around junction anchor)
    # =========================================================
    if ejc:
        _, jpos, _ = ejc
        out["ejc_pm100"] = _extract_seq(genome, chrom, jpos-100, jpos+100, strand)

    return out
    


# ---------- multiprocessing-safe worker setup ----------
_GENOME = None
_DB = None
def _init_worker(genome_path: str, db_path: str):
    global _GENOME, _DB
    _GENOME = Fasta(genome_path)
    _DB = gffutils.FeatureDB(db_path)

def _variant_worker(args):
    idx, chrom, pos, txnames = args
    try:
        seqs = _build_regions_for_variant(_GENOME, _DB, idx, chrom, pos, txnames)
        return idx, seqs, []
    except Exception as e:
        return idx, {}, [str(e)]


# ---------- FIMO runner ----------
def _run_fimo_batch(args):
    batch_id, seq_items, motif_file, threshold, out_base, timeout, is_rna = args
    bdir = out_base / f"batch_{batch_id:05d}"
    rdir = bdir / "results"
    bdir.mkdir(parents=True, exist_ok=True)

    fasta = bdir / "seqs.fasta"
    with open(fasta, "w") as fh:
        for key, seq in seq_items:
            if not seq: continue
            if is_rna:
                seq = seq.replace('T','U').replace('t','u')
            fh.write(f">{key}\n{seq}\n")
    cmd = ["fimo", "--oc", str(rdir), "--thresh", str(threshold), str(motif_file), str(fasta)]
    try:
        subprocess.run(cmd, check=True, timeout=timeout, capture_output=True, text=True)
        tsv = rdir / "fimo.tsv"
        if tsv.exists():
            df = pd.read_csv(tsv, sep="\t", comment="#")
        else:
            df = pd.DataFrame()
        return batch_id, df, ""
    except subprocess.TimeoutExpired:
        return batch_id, pd.DataFrame(), f"timeout_{timeout}s"
    except subprocess.CalledProcessError as e:
        return batch_id, pd.DataFrame(), e.stderr or str(e)

# ---------- Pipeline ----------
class PTCPipelineRFaithful:
    REGION_ORDER = ["utr3_all", "utr3_200","ptc_to_ejc", "newutr_all" , "newutr_200", "ptc_pm100","ejc_pm100"]
    def __init__(self, working_dir=WORKDIR):
        print("="*70); print("PTC REGIONS MOTIF PIPELINE (R-faithful)"); print("="*70)
        self.workdir = Path(working_dir)
        self.variant_file = VARIANTS
        self.gtf_file     = GTF
        self.genome_file  = GENOME
        self.motif_file   = MOTIFS
        self.attract_file = ATTRACT

        for name, path in [
            ("Variants CSV", self.variant_file),
            ("GTF", self.gtf_file),
            ("Genome", self.genome_file),
            ("Motifs", self.motif_file),
            ("ATtRACT", self.attract_file),
        ]:
            ok = Path(path).exists()
            print(f"{name}: {'✓' if ok else '✗ MISSING'} -> {path}")
            if not ok and name != "ATtRACT":
                raise FileNotFoundError(f"{name} missing at {path}")

        print("\nChecking FIMO ..."); print(" ", _check_fimo())
        self.is_rna = IS_RNA

    def setup(self):
        print("\n" + "="*70); print("SETUP"); print("="*70)
        print("Loading genome...")
        self.genome = Fasta(str(self.genome_file))
        print(f"Genome loaded: {len(list(self.genome.keys()))} sequences")

        self.db_path = str(self.gtf_file) + ".txdb.sqlite"
        if Path(self.db_path).exists():
            print("Loading gffutils DB...")
            self.db = gffutils.FeatureDB(self.db_path)
            print("DB loaded.")
        else:
            print("Building gffutils DB (one-time)...")
            gffutils.create_db(
                str(self.gtf_file), dbfn=self.db_path, force=True,
                disable_infer_transcripts=True, disable_infer_genes=True, merge_strategy="merge"
            )
            self.db = gffutils.FeatureDB(self.db_path)
            print("DB built.")

    def load_variants(self, test_mode=False, n_test=2000):
        print("\n" + "="*70); print("LOADING VARIANTS"); print("="*70)
        df = pd.read_csv(self.variant_file)
        # Preserve row order and columns — DO NOT modify the frame beyond filtering stop_gained if present
        if "Consequence" in df.columns:
            mask = df["Consequence"].astype(str).str.contains("stop_gained", case=False, regex=False)
            df = df.loc[mask].copy()
        df = df.reset_index(drop=True)
        if test_mode and len(df) > n_test:
            df = df.iloc[:n_test].copy()
        print(f"Rows: {len(df):,}")
        self.df = df

    def build_regions(self):
        print("\n" + "="*70); print("BUILDING REGIONS & SEQUENCES"); print("="*70)
        args = [(i, self.df.loc[i,"CHROM"], int(self.df.loc[i,"position"]), self.df.loc[i,"txnames"])
                for i in range(len(self.df))]
        self.seq_by_region = {r:{} for r in self.REGION_ORDER}
        errors = []
        for idx, chrom, pos, txnames in args:
            seqdict = _build_regions_for_variant(self.genome, self.db, idx, chrom, pos, txnames)
            print("Returned regions:", seqdict.keys())
            for k, v in seqdict.items():
                if k in self.seq_by_region and v:
                    self.seq_by_region[k][idx] = v
                    errors.extend(errors)
            if errors:
                print(f"Region build warnings: {len(errors)} (showing first 5):", errors[:5])
    
    
    def run_fimo(self):
        print("\n" + "="*70); print("RUNNING FIMO"); print("="*70)
        self.fimo_by_region = {}
        base_out = self.workdir / "TOPMed_Feb_fimo_regions"
        base_out.mkdir(exist_ok=True)
        for region in self.REGION_ORDER:
            seqs = self.seq_by_region.get(region, {})
            if not seqs:
                print(f"[{region}] skipped (no sequences).")
                self.fimo_by_region[region] = pd.DataFrame(); continue
            items = sorted(seqs.items())                  # (row_idx, seq)
            items = [(f"seq_{i}", s) for i, s in items]   # stable header mapping
            n_batches = (len(items) + FIMO_BATCH_SIZE - 1)//FIMO_BATCH_SIZE
            print(f"[{region}] sequences={len(items):,}  batches={n_batches}")
            out_dir = base_out / region; out_dir.mkdir(exist_ok=True)
            batches = []
            for b in range(n_batches):
                start = b*FIMO_BATCH_SIZE; end = min((b+1)*FIMO_BATCH_SIZE, len(items))
                batches.append((b, items[start:end], self.motif_file, FIMO_THRESHOLD, out_dir, FIMO_TIMEOUT_S, self.is_rna))
            with Pool(min(N_CORES, n_batches)) as pool:
                outs = pool.map(_run_fimo_batch, batches)
            dfs, errs = [], 0
            for _, df, err in outs:
                if err: errs += 1
                if not df.empty: dfs.append(df)
            if errs: print(f"[{region}] batches with issues: {errs}")
            self.fimo_by_region[region] = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
            print(f"[{region}] hits: {len(self.fimo_by_region[region]):,}")

    def build_matrices(self):
        print("\n" + "="*70); print("BUILDING MOTIF MATRICES (ATtRACT mapping)"); print("="*70)
        attract = pd.read_csv(self.attract_file, sep="\t")
        # Primary mapping: exact Matrix_id -> Gene_name
        id_to_gene = {}
        if set(["Matrix_id","Gene_name"]).issubset(set(attract.columns)):
            for _, row in attract[["Matrix_id","Gene_name"]].dropna().iterrows():
                id_to_gene[str(row["Matrix_id"])] = str(row["Gene_name"])
        genes_from_ids = set(id_to_gene.values())

        # Fallback: parse gene from motif_id strings like "<id>:GENE:..."
        gene_from_motifid = {}
        pattern = re.compile(r"^[^:]+:([A-Za-z0-9_]+):")
        # gather all motif_ids to seed possible genes
        all_motif_ids = set()
        for region, df in self.fimo_by_region.items():
            if not df.empty and "motif_id" in df:
                all_motif_ids.update(df["motif_id"].astype(str).unique().tolist())
        for mid in all_motif_ids:
            m = pattern.match(mid)
            if m:
                gene_from_motifid[mid] = m.group(1)

        genes = sorted(genes_from_ids.union(set(gene_from_motifid.values())))
        n = len(self.df)

        matrices = []
        for region in self.REGION_ORDER:
            fimo_df = self.fimo_by_region.get(region, pd.DataFrame())
            mat = pd.DataFrame(0, index=np.arange(n), columns=[f"{region}.{g}" for g in genes], dtype=np.int8)
            if not fimo_df.empty and "sequence_name" in fimo_df:
                seq_idx = fimo_df["sequence_name"].astype(str).str.replace("^seq_","", regex=True)
                seq_idx = pd.to_numeric(seq_idx, errors="coerce").astype("Int64")
                fimo_df = fimo_df.assign(seq_idx=seq_idx).dropna(subset=["seq_idx"])

                # Map motif -> gene
                def _to_gene(mid):
                    gid = id_to_gene.get(str(mid))
                    if gid: return gid
                    # try prefix before colon
                    pre = str(mid).split(":")[0]
                    gid = id_to_gene.get(pre)
                    if gid: return gid
                    # parse fallback
                    return gene_from_motifid.get(str(mid), "")

                fimo_df["Gene_name"] = fimo_df["motif_id"].apply(_to_gene).fillna("")
                for g in genes:
                    hits = fimo_df.loc[fimo_df["Gene_name"] == g, "seq_idx"].dropna().astype(int).values
                    if hits.size:
                        hits = np.unique(hits[(hits >= 0) & (hits < n)])
                        if hits.size:
                            mat.loc[hits, f"{region}.{g}"] = 1
            matrices.append(mat)

        full = pd.concat(matrices, axis=1)
        outdir = self.workdir / "topmed_Feb_matrices"; outdir.mkdir(exist_ok=True)
        full.to_csv(outdir / "topmed_motif_presence.csv", index_label="variant_idx")
        print("Saved matrix ->", outdir / "topmed_motif_presence.csv")
        self.matrix = full

    def run(self, test_mode=False):
        self.setup()
        self.load_variants(test_mode=test_mode)
        self.build_regions()
        self.run_fimo()
        self.build_matrices()
        return self.matrix

if __name__ == "__main__":
    pipe = PTCPipelineRFaithful(WORKDIR)
    # Start with a small test (like your R code loops did); flip to False for full.
    pipe.run(test_mode=False)
