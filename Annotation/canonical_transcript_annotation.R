# ==============================================================================
# Variant annotation and canonical transcript selection
#
# Purpose:
#   1. Read ANNOVAR exonic variant annotations
#   2. Construct standardized genomic variant identifiers
#   3. Retain stop-gain SNVs
#   4. Merge ANNOVAR annotations with the dataset-specific variant table
#   5. identify the canonical transcript used for downstream feature generation
#
# Input:
#   - ANNOVAR *.exonic_variant_function file
#   - standardized variant table (`fr.var.can`)
#   - canonical transcript reference table
#
# Output:
#   - variant.anno.hiqual
#     Variants with matched canonical transcript annotations
#
# Genome build:
#   GRCh38 / hg38
#
# Transcript annotation:
#   ANNOVAR ensGene
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)


# ------------------------------------------------------------------------------
# 2. Configuration
# ------------------------------------------------------------------------------

CONFIG <- list(

    annovar_file =
        "/path/to/TOPMed_PTC_SNV_ensGene.exonic_variant_function",

    variant_file =
        "/path/to/fr.var.can.RData",

    canonical_transcript_file =
        "/path/to/canonical_transcripts.tsv",

    output_file =
        "/path/to/TOPMed_canonical_variant_annotation.rds"
)


# ------------------------------------------------------------------------------
# 3. Read ANNOVAR exonic annotations
# ------------------------------------------------------------------------------

variant_anno <- fread(
    CONFIG$annovar_file,
    data.table = FALSE
)

message(
    "ANNOVAR records loaded: ",
    nrow(variant_anno)
)


# ------------------------------------------------------------------------------
# 4. Construct standardized variant key
# ------------------------------------------------------------------------------

variant_anno$key <- paste0(
    variant_anno$V4,
    ":",
    variant_anno$V5,
    "_",
    variant_anno$V7,
    ">",
    variant_anno$V8
)


# ------------------------------------------------------------------------------
# 5. Retain columns required downstream
# ------------------------------------------------------------------------------

variant_anno <- variant_anno[
    ,
    c(
        "V2",
        "V3",
        "V4",
        "V5",
        "V6",
        "V7",
        "V8",
        "key"
    )
]

variant_anno <- unique(
    variant_anno
)


# ------------------------------------------------------------------------------
# 6. Retain stop-gain variants
# ------------------------------------------------------------------------------

variant_anno <- variant_anno[
    grepl(
        "stopgain",
        variant_anno$V2,
        ignore.case = TRUE
    ),
]

message(
    "Stop-gain annotation records retained: ",
    nrow(variant_anno)
)


# ------------------------------------------------------------------------------
# 7. Load standardized variant/genotype table
# ------------------------------------------------------------------------------

load(
    CONFIG$variant_file
)

# Expected object:
#   fr.var.can
#
# Recommended future improvement:
# replace load() with readRDS() so that the object name is explicit.


if (!exists("fr.var.can")) {
    stop(
        "Expected object `fr.var.can` was not found in ",
        CONFIG$variant_file
    )
}


# ------------------------------------------------------------------------------
# 8. Merge ANNOVAR annotation with variant data
# ------------------------------------------------------------------------------

variant_anno_merged <- inner_join(
    variant_anno,
    fr.var.can,
    by = "key"
)

message(
    "Variants matched to genotype/ASE dataset: ",
    nrow(variant_anno_merged)
)


# ------------------------------------------------------------------------------
# 9. Load canonical transcript reference
# ------------------------------------------------------------------------------

canonical_tx <- read.delim(
    CONFIG$canonical_transcript_file,
    stringsAsFactors = FALSE
)

# Expected columns:
#
#   hgnc_symbol
#   ensembl_transcript_id
#
# The reference file should correspond to the transcript annotation
# version used in the manuscript.


# ------------------------------------------------------------------------------
# 10. Extract canonical transcript from ANNOVAR annotation
# ------------------------------------------------------------------------------

get_canonical_transcript <- function(
    annotation,
    canonical_table
) {

    fields <- strsplit(
        annotation,
        ":"
    )[[1]]

    gene <- fields[1]

    canonical_tx <-
        canonical_table %>%
        filter(
            hgnc_symbol == gene
        ) %>%
        pull(
            ensembl_transcript_id
        )

    if (length(canonical_tx) == 0) {
        return(NA_character_)
    }

    # Remove transcript version suffix from ANNOVAR fields for matching
    fields_no_version <- sub(
        "\\..*$",
        "",
        fields
    )

    hits <- which(
        fields_no_version %in%
            canonical_tx
    )

    if (length(hits) != 1) {
        return(NA_character_)
    }

    fields[hits]
}


variant_anno_merged$txnames <- vapply(
    variant_anno_merged$V3,
    get_canonical_transcript,
    character(1),
    canonical_table = canonical_tx
)


# ------------------------------------------------------------------------------
# 11. Retain variants with canonical transcript annotation
# ------------------------------------------------------------------------------

variant.anno.hiqual <-
    variant_anno_merged %>%
    filter(
        !is.na(txnames)
    )


message(
    "Variants with canonical transcript annotation: ",
    nrow(variant.anno.hiqual)
)


# ------------------------------------------------------------------------------
# 12. Save output
# ------------------------------------------------------------------------------

saveRDS(
    variant.anno.hiqual,
    CONFIG$output_file
)

message(
    "Annotation complete: ",
    CONFIG$output_file
)
