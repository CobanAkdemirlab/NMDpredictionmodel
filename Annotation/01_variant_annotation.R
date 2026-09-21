# ==============================================================================
# Shared variant annotation and canonical transcript selection
#
# Purpose:
#   1. Read ANNOVAR exonic variant annotations
#   2. Construct a standardized genomic variant key
#   3. Retain stop-gain annotations
#   4. Merge ANNOVAR annotations with the corresponding standardized
#      dataset-specific variant table
#   5. Select the canonical transcript used for downstream feature generation
#
# Shared across:
#   - TOPMed
#   - gnomAD
#   - ClinVar
#   - GREGoR
#
# Prerequisites:
#   Each dataset must first be processed through Dataset_extraction/ and
#   ANNOVAR. The standardized dataset-specific RDS must contain a `key`
#   column in the format:
#
#       chr7:127588544_A>T
#
# Input:
#   - ANNOVAR *.exonic_variant_function file
#   - standardized dataset-specific PTC SNV RDS
#   - fixed canonical transcript reference table
#
# Output:
#   - <dataset>_canonical_variant_annotation.rds
#
# Genome build:
#   GRCh38 / hg38
#
# ANNOVAR transcript annotation:
#   ensGene
#
# Next step:
#   02_GENCODEv26_transcript_structure.R
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

PROJECT_DIR <- "/path/to/NMDpredictionmodel"

CANONICAL_TRANSCRIPT_FILE <- file.path(
    PROJECT_DIR,
    "Annotation",
    "reference",
    "canonical_transcripts.tsv"
)

OUTPUT_DIR <- file.path(
    PROJECT_DIR,
    "Annotation",
    "output"
)

dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# Dataset-specific inputs.
#
# Each entry uses the SAME downstream processing code.
# Only the paths differ.

DATASETS <- list(

    TOPMed = list(

        annovar_file = file.path(
            PROJECT_DIR,
            "data",
            "TOPMed",
            "TOPMed_PTC_SNV_ensGene.exonic_variant_function"
        ),

        variant_file = file.path(
            PROJECT_DIR,
            "data",
            "TOPMed",
            "TOPMed_PTC_SNV.rds"
        )
    ),


    gnomAD = list(

        annovar_file = file.path(
            PROJECT_DIR,
            "data",
            "gnomAD",
            "gnomAD_v4.1_PTC_SNV_ensGene.exonic_variant_function"
        ),

        variant_file = file.path(
            PROJECT_DIR,
            "data",
            "gnomAD",
            "gnomAD_v4.1_PTC_SNV.rds"
        )
    ),


    ClinVar = list(

        annovar_file = file.path(
            PROJECT_DIR,
            "data",
            "ClinVar",
            "clinvar_20260201_PTC_SNV_ensGene.exonic_variant_function"
        ),

        variant_file = file.path(
            PROJECT_DIR,
            "data",
            "ClinVar",
            "clinvar_20260201_PTC_SNV.rds"
        )
    ),


    GREGoR = list(

        annovar_file = file.path(
            PROJECT_DIR,
            "data",
            "GREGoR",
            "GREGoR_PTC_SNV_ensGene.exonic_variant_function"
        ),

        variant_file = file.path(
            PROJECT_DIR,
            "data",
            "GREGoR",
            "GREGoR_PTC_SNV.rds"
        )
    )
)


# ------------------------------------------------------------------------------
# 3. Load fixed canonical transcript reference
# ------------------------------------------------------------------------------

canonical_tx <- read.delim(
    CANONICAL_TRANSCRIPT_FILE,
    stringsAsFactors = FALSE
)


required_canonical_columns <- c(
    "hgnc_symbol",
    "ensembl_transcript_id"
)

missing_canonical_columns <- setdiff(
    required_canonical_columns,
    colnames(canonical_tx)
)

if (length(missing_canonical_columns) > 0) {

    stop(
        "Canonical transcript file is missing columns: ",
        paste(
            missing_canonical_columns,
            collapse = ", "
        )
    )
}


# Remove transcript version suffix for matching.
#
# Example:
#   ENST00000369636.6
# becomes:
#   ENST00000369636

canonical_tx$ensembl_transcript_id_no_version <- sub(
    "\\..*$",
    "",
    canonical_tx$ensembl_transcript_id
)


# ------------------------------------------------------------------------------
# 4. Canonical transcript selection function
# ------------------------------------------------------------------------------

get_canonical_transcript <- function(
    annotation,
    canonical_table
) {

    if (
        is.na(annotation) ||
        annotation == ""
    ) {
        return(NA_character_)
    }

    # ANNOVAR can contain multiple transcript annotations separated
    # by commas, with individual fields separated by colons.
    fields <- unlist(
        strsplit(
            annotation,
            ":|,"
        )
    )

    if (length(fields) == 0) {
        return(NA_character_)
    }

    # First field is the gene symbol in the ANNOVAR annotation.
    gene <- fields[1]


    # Canonical transcript(s) assigned to this gene.
    canonical_ids <- canonical_table %>%
        filter(
            hgnc_symbol == gene
        ) %>%
        pull(
            ensembl_transcript_id_no_version
        ) %>%
        unique()


    if (length(canonical_ids) == 0) {
        return(NA_character_)
    }


    # Identify ENSEMBL transcript tokens in the ANNOVAR annotation.
    transcript_fields <- fields[
        grepl(
            "^ENST",
            fields
        )
    ]


    if (length(transcript_fields) == 0) {
        return(NA_character_)
    }


    # Remove version suffix only for matching.
    transcript_fields_no_version <- sub(
        "\\..*$",
        "",
        transcript_fields
    )


    hits <- which(
        transcript_fields_no_version %in%
            canonical_ids
    )


    # Require one unambiguous canonical transcript match.
    if (length(hits) != 1) {
        return(NA_character_)
    }


    # Return the transcript exactly as represented by ANNOVAR.
    #
    # Keeping the version suffix, when available, allows direct matching
    # to GENCODE v26 transcript identifiers in subsequent steps.
    transcript_fields[
        hits
    ]
}


# ------------------------------------------------------------------------------
# 5. Process one dataset
# ------------------------------------------------------------------------------

process_dataset <- function(
    dataset_name,
    config
) {

    message(
        "\n============================================================"
    )

    message(
        "Processing dataset: ",
        dataset_name
    )

    message(
        "============================================================"
    )


    # --------------------------------------------------------------------------
    # 5a. Check input files
    # --------------------------------------------------------------------------

    if (!file.exists(config$annovar_file)) {

        stop(
            "ANNOVAR file not found for ",
            dataset_name,
            ": ",
            config$annovar_file
        )
    }


    if (!file.exists(config$variant_file)) {

        stop(
            "Variant file not found for ",
            dataset_name,
            ": ",
            config$variant_file
        )
    }


    # --------------------------------------------------------------------------
    # 5b. Read ANNOVAR exonic annotation
    # --------------------------------------------------------------------------

    variant_anno <- fread(
        config$annovar_file,
        data.table = FALSE
    )


    message(
        "ANNOVAR records loaded: ",
        nrow(variant_anno)
    )


    # ANNOVAR exonic_variant_function columns used here:
    #
    # V2 = functional consequence
    # V3 = gene/transcript/cDNA/protein annotation
    # V4 = chromosome
    # V5 = start
    # V6 = end
    # V7 = reference allele
    # V8 = alternative allele


    # --------------------------------------------------------------------------
    # 5c. Standardize chromosome names
    # --------------------------------------------------------------------------

    variant_anno$V4 <- ifelse(
        grepl(
            "^chr",
            variant_anno$V4
        ),
        variant_anno$V4,
        paste0(
            "chr",
            variant_anno$V4
        )
    )


    # --------------------------------------------------------------------------
    # 5d. Construct standardized variant key
    # --------------------------------------------------------------------------

    variant_anno$key <- paste0(
        variant_anno$V4,
        ":",
        variant_anno$V5,
        "_",
        variant_anno$V7,
        ">",
        variant_anno$V8
    )


    # --------------------------------------------------------------------------
    # 5e. Retain required ANNOVAR columns
    # --------------------------------------------------------------------------

    variant_anno <- variant_anno %>%
        select(
            V2,
            V3,
            V4,
            V5,
            V6,
            V7,
            V8,
            key
        ) %>%
        distinct()


    # --------------------------------------------------------------------------
    # 5f. Retain stop-gain annotations
    # --------------------------------------------------------------------------

    variant_anno <- variant_anno %>%
        filter(
            grepl(
                "stopgain",
                V2,
                ignore.case = TRUE
            )
        )


    message(
        "Stop-gain annotation records retained: ",
        nrow(variant_anno)
    )


    # --------------------------------------------------------------------------
    # 5g. Load standardized variant table
    # --------------------------------------------------------------------------

    variants <- readRDS(
        config$variant_file
    )


    if (!"key" %in% colnames(variants)) {

        stop(
            dataset_name,
            " variant file does not contain a `key` column."
        )
    }


    message(
        "Dataset-specific PTC SNVs loaded: ",
        nrow(variants)
    )


    # --------------------------------------------------------------------------
    # 5h. Merge ANNOVAR annotations with source variant information
    # --------------------------------------------------------------------------

    variant_anno_merged <- inner_join(
        variant_anno,
        variants,
        by = "key"
    )


    message(
        "Records matched between ANNOVAR and source dataset: ",
        nrow(variant_anno_merged)
    )


    # --------------------------------------------------------------------------
    # 5i. Select canonical transcript
    # --------------------------------------------------------------------------

    variant_anno_merged$txnames <- vapply(
        variant_anno_merged$V3,
        get_canonical_transcript,
        character(1),
        canonical_table = canonical_tx
    )


    # --------------------------------------------------------------------------
    # 5j. Retain variants with canonical transcript annotation
    # --------------------------------------------------------------------------

    variants_canonical <- variant_anno_merged %>%
        filter(
            !is.na(txnames)
        )


    message(
        "Records with canonical transcript: ",
        nrow(variants_canonical)
    )


    # --------------------------------------------------------------------------
    # 5k. Add dataset identifier
    # --------------------------------------------------------------------------

    variants_canonical$dataset <- dataset_name


    # --------------------------------------------------------------------------
    # 5l. Save dataset-specific output
    # --------------------------------------------------------------------------

    output_file <- file.path(
        OUTPUT_DIR,
        paste0(
            dataset_name,
            "_canonical_variant_annotation.rds"
        )
    )


    saveRDS(
        variants_canonical,
        output_file
    )


    message(
        "Saved: ",
        output_file
    )


    # --------------------------------------------------------------------------
    # 5m. QC summary
    # --------------------------------------------------------------------------

    qc <- data.frame(

        dataset =
            dataset_name,

        source_variants =
            nrow(variants),

        annovar_stopgain_records =
            nrow(variant_anno),

        matched_records =
            nrow(variant_anno_merged),

        canonical_records =
            nrow(variants_canonical),

        stringsAsFactors = FALSE
    )


    return(qc)
}


# ------------------------------------------------------------------------------
# 6. Run shared workflow across all datasets
# ------------------------------------------------------------------------------

qc_results <- lapply(
    names(DATASETS),
    function(dataset_name) {

        process_dataset(
            dataset_name =
                dataset_name,

            config =
                DATASETS[
                    [dataset_name]
                ]
        )
    }
)


qc_results <- bind_rows(
    qc_results
)


# ------------------------------------------------------------------------------
# 7. Save QC summary
# ------------------------------------------------------------------------------

qc_file <- file.path(
    OUTPUT_DIR,
    "variant_annotation_QC_summary.tsv"
)


write.table(
    qc_results,
    qc_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


print(
    qc_results
)


message(
    "\nShared variant annotation complete."
)

message(
    "Next step:"
)

message(
    "  02_GENCODEv26_transcript_structure.R"
)
