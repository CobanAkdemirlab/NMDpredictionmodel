"""
transcript_to_gene_mapping.py

Purpose
-------
Create a transcript-to-gene mapping from a GENCODE v26 GTF file and
optionally restrict the mapping to transcripts present in a variant dataset.

This helper is used by downstream feature-generation modules that require
conversion between Ensembl transcript IDs (ENST...) and Ensembl gene IDs
(ENSG...).

Reference
---------
Transcript annotation:
    GENCODE v26 primary assembly

Input
-----
1. A CSV containing a transcript identifier column, typically `txnames`
2. GENCODE v26 primary-assembly GTF (.gtf or .gtf.gz)

Output
------
A two-column CSV containing:

    transcript_id
    gene_id

Transcript version suffixes are ignored during matching, while the original
GENCODE transcript and gene identifiers are retained in the output.

This helper is shared across TOPMed, gnomAD, ClinVar, and GREGoR.
"""
import gzip
import re
from pathlib import Path

import pandas as pd


# ==============================================================================
# 1. Configuration
# ==============================================================================

WORKDIR = Path("/path/to/project")

INPUT_CSV = (
    WORKDIR
    / "annotated_variants.csv"
)

GTF_FILE = (
    WORKDIR
    / "reference"
    / "gencode.v26.primary_assembly.annotation.gtf.gz"
)

OUTPUT_FILE = (
    WORKDIR
    / "reference"
    / "transcript_to_gene_mapping_gencode_v26.csv"
)

TRANSCRIPT_COLUMN = "txnames"


# ==============================================================================
# 2. Helper functions
# ==============================================================================

def strip_version(identifier):
    """
    Remove an Ensembl version suffix.

    Example
    -------
    ENST00000368300.8 -> ENST00000368300
    """

    if pd.isna(identifier):
        return None

    return re.sub(
        r"\.\d+$",
        "",
        str(identifier)
    )


def open_gtf(gtf_path):
    """
    Open plain-text or gzip-compressed GTF files.
    """

    if str(gtf_path).endswith(".gz"):
        return gzip.open(
            gtf_path,
            "rt"
        )

    return open(
        gtf_path,
        "rt"
    )


def parse_gencode_tx_gene_mapping(
    gtf_path
):
    """
    Extract transcript_id -> gene_id relationships from a GENCODE GTF.
    """

    mapping = []

    with open_gtf(
        gtf_path
    ) as gtf:

        for line in gtf:

            if line.startswith("#"):
                continue

            fields = (
                line.rstrip("\n")
                .split("\t")
            )

            if len(fields) < 9:
                continue

            attributes = fields[8]

            transcript_match = re.search(
                r'transcript_id "([^"]+)"',
                attributes
            )

            gene_match = re.search(
                r'gene_id "([^"]+)"',
                attributes
            )

            if (
                transcript_match is None
                or
                gene_match is None
            ):
                continue

            transcript_id = (
                transcript_match
                .group(1)
            )

            gene_id = (
                gene_match
                .group(1)
            )

            mapping.append(
                (
                    transcript_id,
                    gene_id
                )
            )

    tx2gene = pd.DataFrame(
        mapping,
        columns=[
            "transcript_id",
            "gene_id"
        ]
    )

    tx2gene = (
        tx2gene
        .drop_duplicates()
        .reset_index(
            drop=True
        )
    )

    return tx2gene


# ==============================================================================
# 3. Main workflow
# ==============================================================================

def main():

    # --------------------------------------------------------------------------
    # Check input files
    # --------------------------------------------------------------------------

    if not INPUT_CSV.exists():

        raise FileNotFoundError(
            f"Input CSV not found: {INPUT_CSV}"
        )

    if not GTF_FILE.exists():

        raise FileNotFoundError(
            f"GENCODE GTF not found: {GTF_FILE}"
        )


    # --------------------------------------------------------------------------
    # Load variant dataset
    # --------------------------------------------------------------------------

    variants = pd.read_csv(
        INPUT_CSV,
        low_memory=False
    )


    if (
        TRANSCRIPT_COLUMN
        not in variants.columns
    ):

        raise ValueError(
            f"Input dataset does not contain "
            f"`{TRANSCRIPT_COLUMN}`."
        )


    transcript_ids = (
        variants[
            TRANSCRIPT_COLUMN
        ]
        .dropna()
        .astype(str)
        .drop_duplicates()
    )


    print(
        "Unique transcript IDs in input:",
        f"{len(transcript_ids):,}"
    )


    # --------------------------------------------------------------------------
    # Parse GENCODE v26 transcript-to-gene mapping
    # --------------------------------------------------------------------------

    tx2gene = (
        parse_gencode_tx_gene_mapping(
            GTF_FILE
        )
    )


    print(
        "Unique transcript-gene pairs in GENCODE v26:",
        f"{len(tx2gene):,}"
    )


    # --------------------------------------------------------------------------
    # Normalize transcript identifiers for matching
    # --------------------------------------------------------------------------

    tx2gene["base_tx_id"] = (
        tx2gene[
            "transcript_id"
        ]
        .map(
            strip_version
        )
    )


    input_transcripts = pd.DataFrame(
        {
            "input_transcript_id":
                transcript_ids
        }
    )


    input_transcripts[
        "base_tx_id"
    ] = (
        input_transcripts[
            "input_transcript_id"
        ]
        .map(
            strip_version
        )
    )


    # --------------------------------------------------------------------------
    # Match input transcripts to GENCODE
    # --------------------------------------------------------------------------

    subset_map = (
        input_transcripts
        .merge(
            tx2gene,
            on="base_tx_id",
            how="left"
        )
    )


    # Count matching based on unique input transcripts rather than rows.
    n_input = (
        input_transcripts[
            "base_tx_id"
        ]
        .nunique()
    )

    n_matched = (
        subset_map[
            "gene_id"
        ]
        .notna()
        .groupby(
            subset_map[
                "base_tx_id"
            ]
        )
        .any()
        .sum()
    )


    print(
        "Matched transcripts:",
        f"{n_matched:,} / {n_input:,}"
    )


    # --------------------------------------------------------------------------
    # Retain one transcript -> gene mapping table
    # --------------------------------------------------------------------------

    output_mapping = (
        subset_map[
            [
                "transcript_id",
                "gene_id"
            ]
        ]
        .dropna()
        .drop_duplicates()
        .sort_values(
            [
                "transcript_id",
                "gene_id"
            ]
        )
        .reset_index(
            drop=True
        )
    )


    # --------------------------------------------------------------------------
    # QC: identify transcripts mapping to multiple genes
    # --------------------------------------------------------------------------

    gene_counts = (
        output_mapping
        .groupby(
            "transcript_id"
        )[
            "gene_id"
        ]
        .nunique()
    )


    ambiguous_transcripts = (
        gene_counts[
            gene_counts > 1
        ]
    )


    if not ambiguous_transcripts.empty:

        print(
            "WARNING:",
            len(
                ambiguous_transcripts
            ),
            "transcripts map to more than one gene."
        )


    # --------------------------------------------------------------------------
    # Save
    # --------------------------------------------------------------------------

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )


    output_mapping.to_csv(
        OUTPUT_FILE,
        index=False
    )


    print(
        "Saved transcript-to-gene mapping:",
        OUTPUT_FILE
    )


# ==============================================================================
# 4. Entry point
# ==============================================================================

if __name__ == "__main__":
    main()
    
