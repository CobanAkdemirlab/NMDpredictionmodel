import os
import math
import pandas as pd
import gffutils
from pyfaidx import Fasta


# Configuration

CSV_FILE = "/NMD_TOPMed/TOPMed_stopgain.csv"
gtf_file = "~/gencode/gencode.v26.primary_assembly.annotation.gtf"
db_file = "~/gencode/gencode_v26.db"

gffutils.create_db(
    gtf_file,
    db_file,
    force=True,
    keep_order=True,
    disable_infer_transcripts=True,
    disable_infer_genes=True
)
GTF_DB = "~/gencode/gencode_v26.db"
print("Database created successfully.")
FASTA_FILE = "~/gencode/hg38.fa"

OUTPUT_DIR = "TOPMed_March_fimo_regions"
BATCH_SIZE = 1000



# Basic Segment Object


class Seg:
    def __init__(self, chrom, start, end, strand):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

    @property
    def width(self):
        return self.end - self.start + 1


# Functions

def parse_cds_exons(cds_exons_str):
    return [int(x) for x in cds_exons_str.split(",") if x.strip()]


def reverse_complement(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def pieces_to_seq(genome, pieces):
    seq = ""
    for chrom, start, end, strand in pieces:
        s = genome[chrom][start-1:end].seq
        if strand == "-":
            s = reverse_complement(s)
        seq += s
    return seq


def order_transcript_direction(features):
    feats = list(features)
    if not feats:
        return []
    strand = feats[0].strand
    feats.sort(key=lambda x: x.start, reverse=(strand == "-"))
    return feats



# CDS coordinate → genomic pieces
def cds_interval_to_genomic(cds_segs, cds_start, cds_end):

    pieces = []
    cds_cursor = 1

    for seg in cds_segs:
        seg_len = seg.width
        seg_cds_start = cds_cursor
        seg_cds_end = cds_cursor + seg_len - 1

        ov_start = max(cds_start, seg_cds_start)
        ov_end = min(cds_end, seg_cds_end)

        if ov_start <= ov_end:

            if seg.strand == "+":
                g_start = seg.start + (ov_start - seg_cds_start)
                g_end = seg.start + (ov_end - seg_cds_start)
            else:
                g_end = seg.end - (ov_start - seg_cds_start)
                g_start = seg.end - (ov_end - seg_cds_start)

            if g_start > g_end:
                g_start, g_end = g_end, g_start

            pieces.append((seg.chrom, int(g_start), int(g_end), seg.strand))

        cds_cursor += seg_len

    return pieces

def slice_firstN_transcript(segs, N):

    result = []
    remaining = N

    for s in segs:
        if remaining <= 0:
            break

        take = min(s.width, remaining)

        if s.strand == "+":
            g_start = s.start
            g_end = s.start + take - 1
        else:
            g_end = s.end
            g_start = s.end - take + 1

        if g_start > g_end:
            g_start, g_end = g_end, g_start

        result.append((s.chrom, g_start, g_end, s.strand))
        remaining -= take

    return result



# Region Builder (CDS-driven)
def build_regions(genome, db, row):

    txid = row["txnames"]
    ptc_cds = int(row["coding.pos"])
    mut_exon = int(row["mut.exon"])
    cds_exons = parse_cds_exons(row["cds_exons"])
    cds_length = int(row["cds_length"])

    ejc_cds = cds_exons[mut_exon - 1]

    try:
        tx = db[txid]
    except:
        return {}

    # CDS segments
    cds_feats = list(db.children(tx, featuretype="CDS", order_by="start"))
    cds_feats = order_transcript_direction(cds_feats)
    cds_segs = [Seg(f.chrom, f.start, f.end, f.strand) for f in cds_feats]

    if not cds_segs:
        return {}

    strand = cds_segs[0].strand

    # 3'UTR segments
    utr_feats = list(db.children(tx, featuretype="UTR", order_by="start"))
    utr_feats = order_transcript_direction(utr_feats)

    # Find CDS "end" in transcript direction
    if strand == "+":
        cds_last = max(f.end for f in cds_feats)   # genomic high end of CDS
        utr_feats = [u for u in utr_feats if u.start > cds_last]
    else:
        cds_last = min(f.start for f in cds_feats) # genomic low end of CDS
        utr_feats = [u for u in utr_feats if u.end < cds_last]

    utr3_segs = [Seg(f.chrom, f.start, f.end, f.strand) for f in utr_feats]

    out = {}

    
    # PTC → EJC
    ptc_to_ejc = cds_interval_to_genomic(cds_segs, ptc_cds, ejc_cds)
    out["ptc_to_ejc"] = pieces_to_seq(genome, ptc_to_ejc)

    
    # NEWUTR_ALL (PTC → CDS end + 3'UTR)
    
    cds_part = cds_interval_to_genomic(cds_segs, ptc_cds, cds_length)
    utr_part = [(s.chrom, s.start, s.end, s.strand) for s in utr3_segs]
    newutr_all = cds_part + utr_part
    out["newutr_all"] = pieces_to_seq(genome, newutr_all)

    # NEWUTR_200
    newutr_200_end = min(ptc_cds + 199, cds_length)
    newutr_200 = cds_interval_to_genomic(cds_segs, ptc_cds, newutr_200_end)
    out["newutr_200"] = pieces_to_seq(genome, newutr_200)

    
    # PTC ±100
   
    ptc_pm = cds_interval_to_genomic(
        cds_segs,
        max(1, ptc_cds - 100),
        min(cds_length, ptc_cds + 100)
    )
    out["ptc_pm100"] = pieces_to_seq(genome, ptc_pm)

    
    # EJC ±100
    
    ejc_pm = cds_interval_to_genomic(
        cds_segs,
        max(1, ejc_cds - 100),
        min(cds_length, ejc_cds + 100)
    )
    out["ejc_pm100"] = pieces_to_seq(genome, ejc_pm)

    
    # UTR3_ALL
    
    utr3_all = [(s.chrom, s.start, s.end, s.strand) for s in utr3_segs]
    out["utr3_all"] = pieces_to_seq(genome, utr3_all)

    # UTR3_200
    utr3_200 = slice_firstN_transcript(utr3_segs, 200)
    out["utr3_200"] = pieces_to_seq(genome, utr3_200)
    #printing out examples 
 

    if row.name in [0, 5]:
        print("\nVariant:", row["key"])
        print("Transcript:", txid)
        print("PTC genomic:", row["position"])
        print("PTC CDS:", ptc_cds)
        print("EJC CDS:", ejc_cds)

        ejc_genomic = cds_interval_to_genomic(cds_segs, ejc_cds, ejc_cds)
        if ejc_genomic:
            chrom, start, end, st = ejc_genomic[0]
            print("EJC genomic:", start, "strand:", st)

        print("strand:", strand)

    return out

# FASTA


def write_fasta_batches(df, region_name):

    region_dir = os.path.join(OUTPUT_DIR, region_name)
    os.makedirs(region_dir, exist_ok=True)

    total = len(df)
    batches = math.ceil(total / BATCH_SIZE)

    for b in range(batches):
        start = b * BATCH_SIZE
        end = min((b + 1) * BATCH_SIZE, total)

        batch_dir = os.path.join(region_dir, f"batch_{b:05d}")
        os.makedirs(batch_dir, exist_ok=True)

        fasta_path = os.path.join(batch_dir, "seqs.fasta")

        with open(fasta_path, "w") as f:
            for i, seq in enumerate(df[region_name].iloc[start:end]):
                if pd.isna(seq) or seq == "":
                    continue
                f.write(f">seq_{start+i}\n")
                f.write(f"{seq}\n")




# MAIN


def main():

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(CSV_FILE, low_memory=False)

    db = gffutils.FeatureDB("~/gencode/gencode_v26.db")
    #remove comment if specific lines needed printing
    #df = df.iloc[[0,5]].copy()

    genome = Fasta(FASTA_FILE)

    results = []

    for _, row in df.iterrows():
        regions = build_regions(genome, db, row)
        results.append(regions)

    regions_df = pd.DataFrame(results)
    final_df = pd.concat([df, regions_df], axis=1)

    # save CSV
    final_df.to_csv("TOPMed_with_sequences.csv", index=False)
    print("Region columns created:")
    print(regions_df.head())
    # aVW FASTA batches for each region
    for region in [
        "ptc_to_ejc",
        "ptc_pm100",
        "ejc_pm100",
        "newutr_all",
        "newutr_200",
        "utr3_all",
        "utr3_200"
    ]:
        write_fasta_batches(final_df, region)
    #printing the examples
    #print("Testing debug rows (original 0 and 5)")

    #for idx in [0, 1]:
        #row = df.iloc[idx]
        #regions = build_regions(genome, db, row)

        #print("\n--- Row positional index:", idx, "original index:", row.name, "---")
        #if not regions:
            #print("No regions returned")
        #else:
            #for k, v in regions.items():
                #print(k, "length =", len(v))

if __name__ == "__main__":
    main()
