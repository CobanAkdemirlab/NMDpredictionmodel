# ==============================================================================
# GENCODE v26 transcript structure extraction
#
# Purpose:
#   Generate transcript-level structural annotations required for downstream
#   PTC positional annotation.
#
# Reference resources:
#   - Genome assembly: GRCh38 / hg38
#   - Transcript annotation: GENCODE v26 primary assembly
#
# Input:
#   - gencode.v26.primary_assembly.annotation.gtf.gz
#
# Output:
#   - gencode_v26_transcript_structure.rds
#
# Output columns:
#   - txnames
#   - cds_length
#   - cds_exons
#   - exon_count
#   - noncoding.length
#   - NCexonsnum
#
# This script is used for TOPMed, gnomAD, ClinVar, and GREGoR.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------------------------

library(GenomicFeatures)
library(AnnotationDbi)
library(GenomicRanges)


# ------------------------------------------------------------------------------
# 2. Configuration
# ------------------------------------------------------------------------------

CONFIG <- list(
    gencode_gtf =
        "/path/to/gencode.v26.primary_assembly.annotation.gtf.gz",

    output_file =
        "/path/to/gencode_v26_transcript_structure.rds"
)


# ------------------------------------------------------------------------------
# 3. Build transcript database from GENCODE v26
# ------------------------------------------------------------------------------

txdb <- makeTxDbFromGFF(
    CONFIG$gencode_gtf
)


# ------------------------------------------------------------------------------
# 4. Restrict to chromosomes used in the analysis
# ------------------------------------------------------------------------------

chr_list <- paste0(
    "chr",
    1:22
)

txdb <- keepSeqlevels(
    txdb,
    chr_list,
    pruning.mode = "coarse"
)


# ------------------------------------------------------------------------------
# 5. Retrieve transcript and CDS annotations
# ------------------------------------------------------------------------------

cds_by_tx <- cdsBy(
    txdb,
    by = "tx",
    use.names = TRUE
)

txnames <- names(
    cds_by_tx
)

tx_metadata <- select(
    txdb,
    keys = txnames,
    columns = columns(txdb),
    keytype = "TXNAME"
)

tx_metadata <- tx_metadata[
    !is.na(tx_metadata$TXNAME),
]


# ------------------------------------------------------------------------------
# 6. CDS length
# ------------------------------------------------------------------------------

get_cds_length <- function(tx) {

    gr <- cds_by_tx[[tx]]

    if (length(gr) == 0) {
        return(NA_real_)
    }

    sum(
        width(gr)
    )
}

cds_length <- vapply(
    txnames,
    get_cds_length,
    numeric(1)
)


# ------------------------------------------------------------------------------
# 7. Cumulative CDS exon boundaries
# ------------------------------------------------------------------------------

get_cds_exons <- function(tx) {

    tx_df <- tx_metadata[
        tx_metadata$TXNAME == tx &
        !is.na(tx_metadata$CDSID),
    ]

    if (nrow(tx_df) == 0) {
        return(NA_character_)
    }

    cds_exon_size <-
        tx_df$CDSEND -
        tx_df$CDSSTART +
        1

    cumulative_exon_ends <-
        cumsum(
            cds_exon_size
        )

    paste(
        cumulative_exon_ends,
        collapse = ","
    )
}

cds_exons <- vapply(
    txnames,
    get_cds_exons,
    character(1)
)


# ------------------------------------------------------------------------------
# 8. Number of coding exons
# ------------------------------------------------------------------------------

get_exon_count <- function(tx) {

    tx_df <- tx_metadata[
        tx_metadata$TXNAME == tx &
        !is.na(tx_metadata$CDSID),
    ]

    nrow(
        tx_df
    )
}

exon_count <- vapply(
    txnames,
    get_exon_count,
    integer(1)
)


# ------------------------------------------------------------------------------
# 9. Total noncoding exon sequence length
# ------------------------------------------------------------------------------

get_noncoding_length <- function(tx) {

    tx_df <- tx_metadata[
        tx_metadata$TXNAME == tx &
        is.na(tx_metadata$CDSID),
    ]

    if (nrow(tx_df) == 0) {
        return(0)
    }

    exon_size <-
        tx_df$EXONEND -
        tx_df$EXONSTART +
        1

    sum(
        exon_size
    )
}

noncoding_length <- vapply(
    txnames,
    get_noncoding_length,
    numeric(1)
)


# ------------------------------------------------------------------------------
# 10. Number of noncoding exons
# ------------------------------------------------------------------------------

get_noncoding_exon_count <- function(tx) {

    tx_df <- tx_metadata[
        tx_metadata$TXNAME == tx &
        is.na(tx_metadata$CDSID),
    ]

    nrow(
        tx_df
    )
}

NCexonsnum <- vapply(
    txnames,
    get_noncoding_exon_count,
    integer(1)
)


# ------------------------------------------------------------------------------
# 11. Build transcript structure table
# ------------------------------------------------------------------------------

transcript_structure <- data.frame(
    txnames = txnames,
    cds_length = cds_length,
    cds_exons = cds_exons,
    exon_count = exon_count,
    noncoding.length = noncoding_length,
    NCexonsnum = NCexonsnum,
    stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# 12. Save output
# ------------------------------------------------------------------------------

saveRDS(
    transcript_structure,
    CONFIG$output_file
)

message(
    "GENCODE v26 transcript structure extraction complete."
)

message(
    "Transcripts retained: ",
    nrow(transcript_structure)
)
