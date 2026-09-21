# ==============================================================================
# Shared variant filtering for downstream analysis
#
# Purpose:
#   Apply the shared biological and expression-based eligibility filters used
#   before downstream feature generation and prediction.
#
# Shared across:
#   - TOPMed
#   - gnomAD
#   - ClinVar
#   - GREGoR
#
# Filtering criteria:
#   1. Retain multi-exon transcripts only
#   2. Remove genes with Whole Blood median expression <= 1 TPM
#   3. Remove variants corresponding to significant Whole Blood eGenes/eQTLs
#   4. Remove variants in known imprinted genes
#
# Important:
#   This script does NOT:
#     - calculate ALLELE.RAT
#     - define NMD.ESCAPEE
#     - randomly select recurrent carriers
#     - perform TOPMed ASE simulation
#
# Those TOPMed-specific steps are handled separately in:
#   05_TOPMed_ASE_simulation.R
#
# Input:
#   - <dataset>_PTC_annotated.rds
#   - GTEx v8 gene median TPM table
#   - GTEx v8 Whole Blood eGenes table
#
# Output:
#   - <dataset>_PTC_filtered.rds
#
# Genome build:
#   GRCh38 / hg38
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------------------------

library(dplyr)


# ------------------------------------------------------------------------------
# 2. Configuration
# ------------------------------------------------------------------------------

CONFIG <- list(

    # Dataset-specific annotated input
    input_file =
        "/path/to/TOPMed_PTC_annotated.rds",

    # GTEx v8 gene median TPM
    gtex_expression_file =
        "/path/to/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.txt",

    # GTEx v8 Whole Blood eGenes
    gtex_egenes_file =
        "/path/to/Whole_Blood.v8.egenes.txt",

    # Output
    output_file =
        "/path/to/TOPMed_PTC_filtered.rds",

    # Minimum Whole Blood median expression
    min_whole_blood_tpm = 1
)


# ------------------------------------------------------------------------------
# 3. Load annotated variants
# ------------------------------------------------------------------------------

variants <- readRDS(
    CONFIG$input_file
)

message(
    "Input annotated variants: ",
    nrow(variants)
)


# ------------------------------------------------------------------------------
# 4. Check required columns
# ------------------------------------------------------------------------------

required_columns <- c(
    "key",
    "variantID",
    "GENE_ID",
    "exon_count"
)

missing_columns <- setdiff(
    required_columns,
    colnames(variants)
)

if (length(missing_columns) > 0) {

    stop(
        "Input dataset is missing required columns: ",
        paste(
            missing_columns,
            collapse = ", "
        )
    )
}


# ------------------------------------------------------------------------------
# 5. Retain multi-exon transcripts only
# ------------------------------------------------------------------------------

n_before <- nrow(variants)

variants <- variants %>%
    filter(
        !is.na(exon_count),
        exon_count > 1
    )

message(
    "Removed single-exon / unresolved transcripts: ",
    n_before - nrow(variants)
)

message(
    "Variants remaining after multi-exon filter: ",
    nrow(variants)
)


# ------------------------------------------------------------------------------
# 6. Remove genes with low Whole Blood expression
#
# Genes with Whole Blood median TPM <= 1 are excluded.
# ------------------------------------------------------------------------------

gtex_expression <- read.table(
    CONFIG$gtex_expression_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# Remove Ensembl version suffix where present
gtex_expression$gene_id <- sub(
    "\\..*$",
    "",
    gtex_expression$Name
)

# Whole Blood column name may be converted by read.table() if check.names=TRUE.
# With check.names=FALSE, retain the original GTEx column name if available.
whole_blood_col <- c(
    "Whole Blood",
    "Whole.Blood"
)

whole_blood_col <- whole_blood_col[
    whole_blood_col %in%
        colnames(gtex_expression)
]

if (length(whole_blood_col) == 0) {

    stop(
        "Could not identify Whole Blood expression column in GTEx table."
    )
}

whole_blood_col <- whole_blood_col[1]

low_expression_genes <- gtex_expression %>%
    filter(
        .data[[whole_blood_col]] <=
            CONFIG$min_whole_blood_tpm
    ) %>%
    pull(
        Description
    ) %>%
    unique()


n_before <- nrow(variants)

variants <- variants %>%
    filter(
        !GENE_ID %in%
            low_expression_genes
    )

message(
    "Removed variants in genes with Whole Blood TPM <= ",
    CONFIG$min_whole_blood_tpm,
    ": ",
    n_before - nrow(variants)
)

message(
    "Variants remaining after Whole Blood expression filter: ",
    nrow(variants)
)


# ------------------------------------------------------------------------------
# 7. Remove significant Whole Blood eQTL-associated variants
#
# GTEx Whole Blood eGenes with q-value <= 0.05 are used.
# The corresponding lead variant IDs are converted to the same variantID
# format used in the analysis by removing the trailing "_b38".
# ------------------------------------------------------------------------------

whole_blood_egenes <- read.table(
    CONFIG$gtex_egenes_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)


required_egene_columns <- c(
    "qval",
    "variant_id"
)

missing_egene_columns <- setdiff(
    required_egene_columns,
    colnames(whole_blood_egenes)
)

if (length(missing_egene_columns) > 0) {

    stop(
        "GTEx Whole Blood eGenes file is missing columns: ",
        paste(
            missing_egene_columns,
            collapse = ", "
        )
    )
}


significant_egenes <- whole_blood_egenes %>%
    filter(
        qval <= 0.05
    )


significant_egenes$variantID <- sub(
    "_b38$",
    "",
    significant_egenes$variant_id
)


eqtl_variant_ids <- unique(
    significant_egenes$variantID
)


n_before <- nrow(variants)

variants <- variants %>%
    filter(
        !variantID %in%
            eqtl_variant_ids
    )

message(
    "Removed significant Whole Blood eQTL-associated variants: ",
    n_before - nrow(variants)
)

message(
    "Variants remaining after eQTL filter: ",
    nrow(variants)
)


# ------------------------------------------------------------------------------
# 8. Remove imprinted genes
# ------------------------------------------------------------------------------

imprinted_genes <- c(
    "UTS2",
    "PPIEL",
    "INPP5F_V2",
    "H19",
    "IGF2",
    "KCNQ1",
    "LPAR6",
    "MEG3",
    "RP11-7F17.7",
    "SNRPN",
    "SNHG14",
    "SNURF",
    "UBE3A",
    "ZNF597",
    "ZNF331",
    "L3MBTL1",
    "NAP1L5",
    "FAM50B",
    "PLAGL1",
    "GRB10",
    "MEST"
)


n_before <- nrow(variants)

variants <- variants %>%
    filter(
        !GENE_ID %in%
            imprinted_genes
    )

message(
    "Removed variants in imprinted genes: ",
    n_before - nrow(variants)
)

message(
    "Variants remaining after imprinting filter: ",
    nrow(variants)
)


# ------------------------------------------------------------------------------
# 9. Final QC summary
# ------------------------------------------------------------------------------

message(
    "\nShared filtering complete."
)

message(
    "Final variants retained: ",
    nrow(variants)
)

message(
    "Unique genomic variants retained: ",
    dplyr::n_distinct(
        variants$key
    )
)

message(
    "Unique genes retained: ",
    dplyr::n_distinct(
        variants$GENE_ID
    )
)


# ------------------------------------------------------------------------------
# 10. Save filtered dataset
# ------------------------------------------------------------------------------

variants_filtered <- variants

saveRDS(
    variants_filtered,
    CONFIG$output_file
)

message(
    "Saved filtered dataset: ",
    CONFIG$output_file
)
