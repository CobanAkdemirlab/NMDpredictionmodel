# ==============================================================================
# GREGoR PTC SNV extraction
#
# Purpose:
#   Identify premature termination codon (PTC)-producing SNVs in the
#   GREGoR variant dataset and prepare them for ANNOVAR annotation.
#
# Input:
#   filtered_oc_base_matches.unique.vcf
#
# Output:
#   1. GREGoR_PTC_SNV.rds
#   2. GREGoR_PTC_SNV.vcf
#   3. GREGoR_PTC_SNV.avinput
#   4. ANNOVAR ensGene annotation files
#
# Genome build:
#   GRCh38 / hg38
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------

CONFIG <- list(
    input_dir   = "~/GREGoR",
    output_dir  = "~/GREGoR/output",
    annovar_dir = "~/annovar",
    annovar_db  = "~/annovar/tempdir",
    vcf_file    = "filtered_oc_base_matches.unique.vcf"
)

dir.create(
    CONFIG$output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 2. Required packages
# ------------------------------------------------------------------------------

library(aenmd)
library(GenomicRanges)
library(GenomeInfoDb)
library(Biostrings)


# ------------------------------------------------------------------------------
# 3. Read GREGoR VCF
# ------------------------------------------------------------------------------

gregor_vcf_file <- file.path(
    CONFIG$input_dir,
    CONFIG$vcf_file
)

gregor_vcf <- aenmd:::parse_vcf_VariantAnnotation(
    gregor_vcf_file
)

gregor_gr <- gregor_vcf$vcf_rng


# ------------------------------------------------------------------------------
# 4. Standardize chromosome naming for aenmd
# ------------------------------------------------------------------------------

seqlevels(gregor_gr) <- sub(
    "^chr",
    "",
    seqlevels(gregor_gr)
)

seqnames(gregor_gr) <- sub(
    "^chr",
    "",
    as.character(seqnames(gregor_gr))
)


# ------------------------------------------------------------------------------
# 5. Process variants
# ------------------------------------------------------------------------------

gregor_processed <- process_variants(
    gregor_gr
)

# Remove variants containing undefined ALT alleles
contains_N <- Biostrings::vcountPattern(
    "N",
    gregor_processed$alt
) > 0

gregor_processed <- gregor_processed[
    !contains_N
]


# ------------------------------------------------------------------------------
# 6. Annotate NMD / premature termination codons
# ------------------------------------------------------------------------------

gregor_nmd <- annotate_nmd(
    gregor_processed,
    rettype = "gr"
)

gregor_df <- as.data.frame(
    gregor_nmd,
    row.names = NULL,
    optional = TRUE
)


# ------------------------------------------------------------------------------
# 7. Retain PTC-producing variants
# ------------------------------------------------------------------------------

gregor_ptc <- gregor_df[
    gregor_df$res_aenmd.is_ptc %in% TRUE,
]

message(
    "PTC-producing variants before SNV filtering: ",
    nrow(gregor_ptc)
)


# ------------------------------------------------------------------------------
# 8. Retain SNVs only
# ------------------------------------------------------------------------------

gregor_ptc_snv <- gregor_ptc[
    gregor_ptc$type == "snv",
]

message(
    "PTC-producing SNVs retained: ",
    nrow(gregor_ptc_snv)
)


# ------------------------------------------------------------------------------
# 9. Create standardized variant identifiers
# ------------------------------------------------------------------------------

gregor_ptc_snv$CHROM <- as.character(
    gregor_ptc_snv$seqnames
)

gregor_ptc_snv$CHROM <- ifelse(
    grepl("^chr", gregor_ptc_snv$CHROM),
    gregor_ptc_snv$CHROM,
    paste0("chr", gregor_ptc_snv$CHROM)
)

gregor_ptc_snv$POS <- gregor_ptc_snv$start

gregor_ptc_snv$key <- with(
    gregor_ptc_snv,
    paste0(
        CHROM, ":",
        POS, "_",
        ref, ">",
        alt
    )
)

gregor_ptc_snv$variantID <- with(
    gregor_ptc_snv,
    paste(
        CHROM,
        POS,
        ref,
        alt,
        sep = "_"
    )
)


# ------------------------------------------------------------------------------
# 10. Save filtered GREGoR PTC SNVs
# ------------------------------------------------------------------------------

saveRDS(
    gregor_ptc_snv,
    file.path(
        CONFIG$output_dir,
        "GREGoR_PTC_SNV.rds"
    )
)


# ------------------------------------------------------------------------------
# 11. Prepare VCF-like input for ANNOVAR
# ------------------------------------------------------------------------------

annovar_vcf <- data.frame(
    CHROM  = gregor_ptc_snv$CHROM,
    POS    = gregor_ptc_snv$POS,
    ID     = ".",
    REF    = as.character(gregor_ptc_snv$ref),
    ALT    = as.character(gregor_ptc_snv$alt),
    QUAL   = ".",
    FILTER = "PASS",
    INFO   = ".",
    stringsAsFactors = FALSE
)

annovar_vcf_file <- file.path(
    CONFIG$output_dir,
    "GREGoR_PTC_SNV.vcf"
)

write.table(
    annovar_vcf,
    file      = annovar_vcf_file,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
)


# ------------------------------------------------------------------------------
# 12. Convert VCF to ANNOVAR input
# ------------------------------------------------------------------------------

convert2annovar <- file.path(
    CONFIG$annovar_dir,
    "convert2annovar.pl"
)

annovar_avinput <- file.path(
    CONFIG$output_dir,
    "GREGoR_PTC_SNV.avinput"
)

convert_command <- sprintf(
    "perl %s -format vcf4 %s > %s",
    shQuote(convert2annovar),
    shQuote(annovar_vcf_file),
    shQuote(annovar_avinput)
)

status <- system(convert_command)

if (status != 0) {
    stop("convert2annovar.pl failed.")
}


# ------------------------------------------------------------------------------
# 13. Run ANNOVAR ensGene annotation
# ------------------------------------------------------------------------------

annotate_variation <- file.path(
    CONFIG$annovar_dir,
    "annotate_variation.pl"
)

annovar_output_prefix <- file.path(
    CONFIG$output_dir,
    "GREGoR_PTC_SNV_ensGene"
)

annotation_command <- sprintf(
    paste(
        "perl %s",
        "-build hg38",
        "-out %s",
        "-dbtype ensGene",
        "%s",
        "%s"
    ),
    shQuote(annotate_variation),
    shQuote(annovar_output_prefix),
    shQuote(annovar_avinput),
    shQuote(CONFIG$annovar_db)
)

status <- system(annotation_command)

if (status != 0) {
    stop("ANNOVAR gene annotation failed.")
}


# ------------------------------------------------------------------------------
# 14. Completion message
# ------------------------------------------------------------------------------

message("GREGoR PTC SNV extraction complete.")
message("Output directory: ", CONFIG$output_dir)
