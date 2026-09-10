# ==============================================================================
# ClinVar PTC extraction
#
# Purpose:
#   Extract premature termination codon (PTC)-producing SNVs from ClinVar
#   and prepare them for downstream ANNOVAR transcript annotation.
#
# Input:
#   ClinVar GRCh38 VCF:
#       clinvar_20260201.vcf.gz
#
# Output:
#   1. clinvar_20260201_PTC_SNV.rds
#   2. clinvar_20260201_PTC_SNV.vcf
#   3. clinvar_20260201_PTC_SNV.avinput
#   4. ANNOVAR ensGene annotation files
#
# Genome build:
#   GRCh38 / hg38
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------

CONFIG <- list(
    input_dir   = "~/ClinVar2026",
    output_dir  = "~/ClinVar2026/output",
    annovar_dir = "~/annovar",
    annovar_db  = "~/annovar/humandb",
    vcf_file    = "clinvar_20260201.vcf.gz"
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
library(S4Vectors)
library(Biostrings)


# ------------------------------------------------------------------------------
# 3. Read ClinVar VCF
# ------------------------------------------------------------------------------

clinvar_vcf <- file.path(
    CONFIG$input_dir,
    CONFIG$vcf_file
)

clinvar <- aenmd:::parse_vcf_VariantAnnotation(clinvar_vcf)

clinvar_gr <- clinvar$vcf_rng


# ------------------------------------------------------------------------------
# 4. Preprocess variants
# ------------------------------------------------------------------------------

clinvar_processed <- process_variants(clinvar_gr)

# Remove variants containing undefined alternative alleles ("N")
contains_N <- Biostrings::vcountPattern(
    "N",
    clinvar_processed$alt
) > 0

clinvar_processed <- clinvar_processed[!contains_N]


# ------------------------------------------------------------------------------
# 5. Annotate predicted premature termination codons
# ------------------------------------------------------------------------------

clinvar_nmd <- annotate_nmd(
    clinvar_processed,
    rettype = "gr"
)

clinvar_df <- as.data.frame(
    clinvar_nmd,
    row.names = NULL,
    optional = TRUE
)


# ------------------------------------------------------------------------------
# 6. Retain PTC-producing SNVs
# ------------------------------------------------------------------------------

clinvar_ptc <- clinvar_df[
    clinvar_df$res_aenmd.is_ptc %in% TRUE,
]

clinvar_ptc_snv <- clinvar_ptc[
    clinvar_ptc$type == "snv",
]
message("PTC variants before SNV filtering: ", nrow(clinvar_ptc))
message("PTC SNVs retained: ", nrow(clinvar_ptc_snv))

# ------------------------------------------------------------------------------
# 7. Create standardized variant identifiers
# ------------------------------------------------------------------------------

clinvar_ptc_snv$CHROM <- as.character(clinvar_ptc_snv$seqnames)
clinvar_ptc_snv$POS   <- clinvar_ptc_snv$start

# Add "chr" only if chromosome names do not already contain it
clinvar_ptc_snv$CHROM <- ifelse(
    grepl("^chr", clinvar_ptc_snv$CHROM),
    clinvar_ptc_snv$CHROM,
    paste0("chr", clinvar_ptc_snv$CHROM)
)

clinvar_ptc_snv$key <- with(
    clinvar_ptc_snv,
    paste0(CHROM, ":", POS, "_", ref, ">", alt)
)

clinvar_ptc_snv$variantID <- with(
    clinvar_ptc_snv,
    paste(CHROM, POS, ref, alt, sep = "_")
)


# ------------------------------------------------------------------------------
# 8. Save the filtered ClinVar dataset
# ------------------------------------------------------------------------------

saveRDS(
    clinvar_ptc_snv,
    file.path(
        CONFIG$output_dir,
        "clinvar_20260201_PTC_SNV.rds"
    )
)


# ------------------------------------------------------------------------------
# 9. Prepare VCF-like input for ANNOVAR
# ------------------------------------------------------------------------------

annovar_vcf <- data.frame(
    CHROM  = clinvar_ptc_snv$CHROM,
    POS    = clinvar_ptc_snv$POS,
    ID     = ".",
    REF    = clinvar_ptc_snv$ref,
    ALT    = clinvar_ptc_snv$alt,
    QUAL   = ".",
    FILTER = ".",
    INFO   = "."
)

annovar_vcf_file <- file.path(
    CONFIG$output_dir,
    "clinvar_20260201_PTC_SNV.vcf"
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
# 10. Convert VCF to ANNOVAR input
# ------------------------------------------------------------------------------

annovar_avinput <- file.path(
    CONFIG$output_dir,
    "clinvar_20260201_PTC_SNV.avinput"
)

convert2annovar <- file.path(
    CONFIG$annovar_dir,
    "convert2annovar.pl"
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
# 11. Annotate variants using ANNOVAR ensGene
# ------------------------------------------------------------------------------

annotate_variation <- file.path(
    CONFIG$annovar_dir,
    "annotate_variation.pl"
)

annovar_output_prefix <- file.path(
    CONFIG$output_dir,
    "clinvar_20260201_PTC_SNV_ensGene"
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
# 12. Report dataset dimensions
# ------------------------------------------------------------------------------

message(
    "ClinVar PTC extraction complete."
)

message(
    "PTC-producing variants: ",
    nrow(clinvar_ptc)
)

message(
    "PTC-producing SNVs retained: ",
    nrow(clinvar_ptc_snv)
)

message(
    "Output directory: ",
    normalizePath(CONFIG$output_dir)
)


