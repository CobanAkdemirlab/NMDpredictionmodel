"""
07a_motif_region_extraction.py

PTC-centered sequence region extraction for RBP motif analysis.

Purpose
-------
Construct transcript-oriented nucleotide sequence regions surrounding each
premature termination codon (PTC) for subsequent motif scanning with FIMO.

This script is designed to be shared across:
    - TOPMed
    - gnomAD
    - ClinVar
    - GREGoR

Reference resources
-------------------
Genome assembly:
    GRCh38 / hg38

Transcript annotation:
    GENCODE v26 primary assembly

Required Python packages
------------------------
pandas
gffutils
pyfaidx

Required input columns
----------------------
The input variant table must contain:

    key
    txnames
    coding.pos
    mut.exon
    cds_exons
    cds_length

Optional columns such as chromosome or genomic position may also be present
but are not required for region construction.

Regions generated
-----------------
ptc_to_ejc
    Sequence from the PTC to the end of the PTC-bearing coding exon.

ptc_pm100
    CDS sequence spanning up to 100 nt upstream and downstream of the PTC.

ejc_pm100
    CDS sequence spanning up to 100 nt around the end of the PTC-bearing
    coding exon.

newutr_all
    Sequence from the PTC through the remaining CDS plus the annotated 3'UTR.

newutr_200
    First 200 nt beginning at the PTC within the remaining CDS.

utr3_all
    Complete annotated 3'UTR sequence.

utr3_200
    First 200 nt of the annotated 3'UTR.

Output
------
1. A CSV containing the original variant table plus the seven extracted
   sequence regions.

2. Batched FASTA files for each region for input to FIMO.

Next step
---------
07b_FIMO_motif_features.py
"""

import math
from pathlib import Path

import pandas as pd
import gffutils
from pyfaidx import Fasta


# ==============================================================================
# 1. Configuration
# ==============================================================================

WORKDIR = Path("/path/to/project")

INPUT_CSV = WORKDIR / "annotated_variants.csv"

GTF_FILE = (
    WORKDIR
    / "reference"
    / "gencode.v26.primary_assembly.annotation.gtf"
)

GTF_DB = (
    WORKDIR
    / "reference"
    / "gencode_v26.db"
)

FASTA_FILE = (
    WORKDIR
    / "reference"
    / "hg38.fa"
)

OUTPUT_DIR = (
    WORKDIR
    / "motif_regions"
)

OUTPUT_CSV = (
    WORKDIR
    / "variants_with_motif_regions.csv"
)

BATCH_SIZE = 1000

REGIONS = [
    "ptc_to_ejc",
    "ptc_pm100",
    "ejc_pm100",
    "newutr_all",
    "newutr_200",
    "utr3_all",
    "utr3_200",
]


# ==============================================================================
# 2. Basic segment object
# ==============================================================================

class Seg:
    """Simple genomic segment container."""

    def __init__(self, chrom, start, end, strand):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

    @property
    def width(self):
        return self.end - self.start + 1


# ==============================================================================
# 3. Helper functions
# ==============================================================================

def parse_cds_exons(cds_exons_str):
    """
    Parse comma-separated cumulative CDS exon-end coordinates.

    Example
    -------
    "120,245,390" -> [120, 245, 390]
    """

    if pd.isna(cds_exons_str):
        return []

    return [
        int(x)
        for x in str(cds_exons_str).split(",")
        if x.strip()
    ]


def reverse_complement(seq):
    """Return reverse-complement DNA sequence."""

    comp = str.maketrans(
        "ACGTacgt",
        "TGCAtgca"
    )

    return seq.translate(comp)[::-1]


def pieces_to_seq(genome, pieces):
    """
    Convert genomic pieces into one transcript-oriented sequence.

    Each piece is:
        (chrom, start, end, strand)
    """

    seq = ""

    for chrom, start, end, strand in pieces:

        fragment = genome[
            chrom
        ][
            start - 1:end
        ].seq

        if strand == "-":
            fragment = reverse_complement(
                fragment
            )

        seq += fragment

    return seq


def order_transcript_direction(features):
    """
    Sort genomic features in transcript 5' -> 3' direction.
    """

    feats = list(features)

    if not feats:
        return []

    strand = feats[0].strand

    feats.sort(
        key=lambda x: x.start,
        reverse=(strand == "-")
    )

    return feats


# ==============================================================================
# 4. CDS-coordinate to genomic-coordinate conversion
# ==============================================================================

def cds_interval_to_genomic(
    cds_segs,
    cds_start,
    cds_end
):
    """
    Convert a CDS-coordinate interval into one or more genomic intervals.

    Parameters
    ----------
    cds_segs
        CDS segments ordered in transcript direction.

    cds_start
        1-based CDS start coordinate.

    cds_end
        1-based CDS end coordinate.

    Returns
    -------
    list of tuples:
        (chrom, genomic_start, genomic_end, strand)
    """

    pieces = []
    cds_cursor = 1

    for seg in cds_segs:

        seg_len = seg.width

        seg_cds_start = cds_cursor
        seg_cds_end = (
            cds_cursor +
            seg_len -
            1
        )

        ov_start = max(
            cds_start,
            seg_cds_start
        )

        ov_end = min(
            cds_end,
            seg_cds_end
        )

        if ov_start <= ov_end:

            if seg.strand == "+":

                g_start = (
                    seg.start +
                    (
                        ov_start -
                        seg_cds_start
                    )
                )

                g_end = (
                    seg.start +
                    (
                        ov_end -
                        seg_cds_start
                    )
                )

            else:

                g_end = (
                    seg.end -
                    (
                        ov_start -
                        seg_cds_start
                    )
                )

                g_start = (
                    seg.end -
                    (
                        ov_end -
                        seg_cds_start
                    )
                )

            if g_start > g_end:
                g_start, g_end = (
                    g_end,
                    g_start
                )

            pieces.append(
                (
                    seg.chrom,
                    int(g_start),
                    int(g_end),
                    seg.strand
                )
            )

        cds_cursor += seg_len

    return pieces


def slice_first_n_transcript(
    segs,
    n
):
    """
    Return genomic pieces corresponding to the first n nucleotides
    of a transcript-oriented segmented region.
    """

    result = []
    remaining = n

    for seg in segs:

        if remaining <= 0:
            break

        take = min(
            seg.width,
            remaining
        )

        if seg.strand == "+":

            g_start = seg.start
            g_end = (
                seg.start +
                take -
                1
            )

        else:

            g_end = seg.end
            g_start = (
                seg.end -
                take +
                1
            )

        if g_start > g_end:
            g_start, g_end = (
                g_end,
                g_start
            )

        result.append(
            (
                seg.chrom,
                g_start,
                g_end,
                seg.strand
            )
        )

        remaining -= take

    return result


# ==============================================================================
# 5. GENCODE database preparation
# ==============================================================================

def prepare_gencode_db():
    """
    Create the gffutils GENCODE v26 database if it does not already exist.
    """

    GTF_DB.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    if GTF_DB.exists():

        print(
            "Using existing GENCODE v26 database:",
            GTF_DB
        )

        return

    print(
        "Creating GENCODE v26 database..."
    )

    gffutils.create_db(
        str(GTF_FILE),
        str(GTF_DB),
        force=False,
        keep_order=True,
        disable_infer_transcripts=True,
        disable_infer_genes=True
    )

    print(
        "GENCODE database created:",
        GTF_DB
    )


# ==============================================================================
# 6. Region builder
# ==============================================================================

def build_regions(
    genome,
    db,
    row
):
    """
    Construct the seven sequence regions used for motif analysis.
    """

    required_fields = [
        "txnames",
        "coding.pos",
        "mut.exon",
        "cds_exons",
        "cds_length",
    ]

    for field in required_fields:

        if field not in row.index:
            raise KeyError(
                f"Required column missing: {field}"
            )

    if (
        pd.isna(row["txnames"]) or
        pd.isna(row["coding.pos"]) or
        pd.isna(row["mut.exon"]) or
        pd.isna(row["cds_exons"]) or
        pd.isna(row["cds_length"])
    ):
        return {}

    txid = str(
        row["txnames"]
    )

    ptc_cds = int(
        row["coding.pos"]
    )

    mut_exon = int(
        row["mut.exon"]
    )

    cds_exons = parse_cds_exons(
        row["cds_exons"]
    )

    cds_length = int(
        row["cds_length"]
    )

    if not cds_exons:
        return {}

    if (
        mut_exon < 1 or
        mut_exon > len(cds_exons)
    ):
        return {}

    # End of PTC-bearing coding exon in CDS coordinates
    ejc_cds = cds_exons[
        mut_exon - 1
    ]

    try:
        tx = db[txid]

    except Exception:
        return {}

    # --------------------------------------------------------------------------
    # CDS segments
    # --------------------------------------------------------------------------

    cds_feats = list(
        db.children(
            tx,
            featuretype="CDS",
            order_by="start"
        )
    )

    cds_feats = order_transcript_direction(
        cds_feats
    )

    cds_segs = [
        Seg(
            f.chrom,
            f.start,
            f.end,
            f.strand
        )
        for f in cds_feats
    ]

    if not cds_segs:
        return {}

    strand = cds_segs[0].strand

    # --------------------------------------------------------------------------
    # 3'UTR segments
    # --------------------------------------------------------------------------

    utr_feats = list(
        db.children(
            tx,
            featuretype="UTR",
            order_by="start"
        )
    )

    utr_feats = order_transcript_direction(
        utr_feats
    )

    # Retain only UTR sequence downstream of the CDS
    # in transcript orientation.

    if strand == "+":

        cds_last = max(
            f.end
            for f in cds_feats
        )

        utr_feats = [
            u
            for u in utr_feats
            if u.start > cds_last
        ]

    else:

        cds_last = min(
            f.start
            for f in cds_feats
        )

        utr_feats = [
            u
            for u in utr_feats
            if u.end < cds_last
        ]

    utr3_segs = [
        Seg(
            f.chrom,
            f.start,
            f.end,
            f.strand
        )
        for f in utr_feats
    ]

    out = {}

    # --------------------------------------------------------------------------
    # Region 1: PTC -> EJC
    # --------------------------------------------------------------------------

    ptc_to_ejc = cds_interval_to_genomic(
        cds_segs,
        ptc_cds,
        ejc_cds
    )

    out["ptc_to_ejc"] = pieces_to_seq(
        genome,
        ptc_to_ejc
    )

    # --------------------------------------------------------------------------
    # Region 2: newUTR_all
    #
    # PTC -> remaining CDS + annotated 3'UTR
    # --------------------------------------------------------------------------

    cds_part = cds_interval_to_genomic(
        cds_segs,
        ptc_cds,
        cds_length
    )

    utr_part = [
        (
            seg.chrom,
            seg.start,
            seg.end,
            seg.strand
        )
        for seg in utr3_segs
    ]

    newutr_all = (
        cds_part +
        utr_part
    )

    out["newutr_all"] = pieces_to_seq(
        genome,
        newutr_all
    )

    # --------------------------------------------------------------------------
    # Region 3: newUTR_200
    #
    # First 200 nt beginning at the PTC within the remaining CDS.
    #
    # NOTE:
    # This preserves the implementation used in the original analysis.
    # It does not extend into the annotated 3'UTR if fewer than 200 CDS
    # nucleotides remain downstream of the PTC.
    # --------------------------------------------------------------------------

    newutr_200_end = min(
        ptc_cds + 199,
        cds_length
    )

    newutr_200 = cds_interval_to_genomic(
        cds_segs,
        ptc_cds,
        newutr_200_end
    )

    out["newutr_200"] = pieces_to_seq(
        genome,
        newutr_200
    )

    # --------------------------------------------------------------------------
    # Region 4: PTC +/-100 nt
    # --------------------------------------------------------------------------

    ptc_pm100 = cds_interval_to_genomic(
        cds_segs,
        max(
            1,
            ptc_cds - 100
        ),
        min(
            cds_length,
            ptc_cds + 100
        )
    )

    out["ptc_pm100"] = pieces_to_seq(
        genome,
        ptc_pm100
    )

    # --------------------------------------------------------------------------
    # Region 5: EJC +/-100 nt
    #
    # Here the EJC-associated coordinate is represented by the end of the
    # PTC-bearing coding exon in CDS coordinates.
    # --------------------------------------------------------------------------

    ejc_pm100 = cds_interval_to_genomic(
        cds_segs,
        max(
            1,
            ejc_cds - 100
        ),
        min(
            cds_length,
            ejc_cds + 100
        )
    )

    out["ejc_pm100"] = pieces_to_seq(
        genome,
        ejc_pm100
    )

    # --------------------------------------------------------------------------
    # Region 6: complete annotated 3'UTR
    # --------------------------------------------------------------------------

    utr3_all = [
        (
            seg.chrom,
            seg.start,
            seg.end,
            seg.strand
        )
        for seg in utr3_segs
    ]

    out["utr3_all"] = pieces_to_seq(
        genome,
        utr3_all
    )

    # --------------------------------------------------------------------------
    # Region 7: first 200 nt of annotated 3'UTR
    # --------------------------------------------------------------------------

    utr3_200 = slice_first_n_transcript(
        utr3_segs,
        200
    )

    out["utr3_200"] = pieces_to_seq(
        genome,
        utr3_200
    )

    return out


# ==============================================================================
# 7. FASTA writer
# ==============================================================================

def make_safe_fasta_id(
    key
):
    """
    Convert standardized variant key into a FASTA-safe identifier.

    Example
    -------
    chr7:127588544_A>T
        ->
    chr7_127588544_A_T
    """

    return (
        str(key)
        .replace(":", "_")
        .replace(">", "_")
        .replace("/", "_")
        .replace(" ", "_")
    )


def write_fasta_batches(
    df,
    region_name
):
    """
    Write batched FASTA files for one sequence region.
    """

    region_dir = (
        OUTPUT_DIR /
        region_name
    )

    region_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    total = len(df)

    batches = math.ceil(
        total /
        BATCH_SIZE
    )

    for batch_index in range(
        batches
    ):

        start = (
            batch_index *
            BATCH_SIZE
        )

        end = min(
            (
                batch_index +
                1
            ) *
            BATCH_SIZE,
            total
        )

        batch_dir = (
            region_dir /
            f"batch_{batch_index:05d}"
        )

        batch_dir.mkdir(
            parents=True,
            exist_ok=True
        )

        fasta_path = (
            batch_dir /
            "seqs.fasta"
        )

        with open(
            fasta_path,
            "w"
        ) as fasta_handle:

            batch_df = df.iloc[
                start:end
            ]

            for _, row in batch_df.iterrows():

                seq = row[
                    region_name
                ]

                if (
                    pd.isna(seq) or
                    seq == ""
                ):
                    continue

                fasta_id = make_safe_fasta_id(
                    row["key"]
                )

                fasta_handle.write(
                    f">{fasta_id}\n"
                )

                fasta_handle.write(
                    f"{seq}\n"
                )


# ==============================================================================
# 8. QC
# ==============================================================================

def validate_input_columns(
    df
):
    """
    Confirm that required columns are available.
    """

    required_columns = [
        "key",
        "txnames",
        "coding.pos",
        "mut.exon",
        "cds_exons",
        "cds_length",
    ]

    missing = [
        col
        for col in required_columns
        if col not in df.columns
    ]

    if missing:

        raise ValueError(
            "Missing required input columns: "
            +
            ", ".join(
                missing
            )
        )


# ==============================================================================
# 9. Main workflow
# ==============================================================================

def main():

    # --------------------------------------------------------------------------
    # Create output directories
    # --------------------------------------------------------------------------

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    # --------------------------------------------------------------------------
    # Prepare GENCODE database
    # --------------------------------------------------------------------------

    prepare_gencode_db()

    # --------------------------------------------------------------------------
    # Load variant dataset
    # --------------------------------------------------------------------------

    print(
        "Loading variants:",
        INPUT_CSV
    )

    variants = pd.read_csv(
        INPUT_CSV,
        low_memory=False
    )

    validate_input_columns(
        variants
    )

    print(
        "Variants loaded:",
        len(variants)
    )

    # --------------------------------------------------------------------------
    # Load GENCODE database and GRCh38 genome
    # --------------------------------------------------------------------------

    db = gffutils.FeatureDB(
        str(GTF_DB)
    )

    genome = Fasta(
        str(FASTA_FILE)
    )

    # --------------------------------------------------------------------------
    # Extract sequence regions
    # --------------------------------------------------------------------------

    results = []

    for index, row in variants.iterrows():

        if (
            index % 500 ==
            0
        ):

            print(
                f"Processing variant "
                f"{index + 1} / "
                f"{len(variants)}"
            )

        regions = build_regions(
            genome,
            db,
            row
        )

        results.append(
            regions
        )

    regions_df = pd.DataFrame(
        results
    )

    final_df = pd.concat(
        [
            variants.reset_index(
                drop=True
            ),
            regions_df.reset_index(
                drop=True
            )
        ],
        axis=1
    )

    # --------------------------------------------------------------------------
    # Save sequence table
    # --------------------------------------------------------------------------

    final_df.to_csv(
        OUTPUT_CSV,
        index=False
    )

    print(
        "Saved sequence-region table:",
        OUTPUT_CSV
    )

    # --------------------------------------------------------------------------
    # QC summary
    # --------------------------------------------------------------------------

    print(
        "\nSequence-region availability:"
    )

    for region in REGIONS:

        if region not in final_df.columns:
            print(
                region,
                ": column missing"
            )

            continue

        n_available = (
            final_df[region]
            .fillna("")
            .ne("")
            .sum()
        )

        print(
            f"  {region}: "
            f"{n_available} / "
            f"{len(final_df)}"
        )

    # --------------------------------------------------------------------------
    # Write FASTA batches
    # --------------------------------------------------------------------------

    print(
        "\nWriting FASTA batches..."
    )

    for region in REGIONS:

        write_fasta_batches(
            final_df,
            region
        )

    print(
        "\nMotif region extraction complete."
    )

    print(
        "Next step:"
    )

    print(
        "  python 07b_FIMO_motif_features.py"
    )


# ==============================================================================
# 10. Entry point
# ==============================================================================

if __name__ == "__main__":
    main()
