# ==============================================================================
# TOPMed ASE/genotype extraction
#
# Purpose:
#   1. Read per-sample TOPMed ASE/genotype files
#   2. Retain heterozygous variants
#   3. Combine variants across individuals
#   4. Construct standardized genomic variant identifiers
#   5. Prepare SNVs for ANNOVAR ensGene annotation
#
# Input:
#   Per-sample TOPMed ASE/genotype VCF-like files generated after
#   ASEReadCounter and genotype matching
#
# Output:
#   TOPMed_heterozygous_PTC_SNVs.tsv
#   TOPMed_PTC_SNV.vcf
#   TOPMed_PTC_SNV.avinput
#   ANNOVAR ensGene annotation files
#
# Genome build:
#   GRCh38 / hg38
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
    input_dir   = "~/ASE_genotype",
    output_dir  = "~/NMD_TOPMed",
    annovar_dir = "~/annovar",
    annovar_db  = "~/annovar/tempdir",
    file_pattern = "\\.vcf$"
)

dir.create(
    CONFIG$output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 3. Locate per-sample files
# ------------------------------------------------------------------------------

files <- list.files(
    CONFIG$input_dir,
    pattern = CONFIG$file_pattern,
    full.names = TRUE
)

if (length(files) == 0) {
    stop("No input VCF files were found in: ", CONFIG$input_dir)
}

message("Number of sample files found: ", length(files))


# ------------------------------------------------------------------------------
# 4. Extract heterozygous variants from each individual
# ------------------------------------------------------------------------------

het_list <- vector(
    "list",
    length(files)
)

for (i in seq_along(files)) {

    message(
        "Processing ",
        i,
        " of ",
        length(files),
        ": ",
        basename(files[i])
    )

    sample_df <- read.table(
        files[i],
        header = TRUE,
        stringsAsFactors = FALSE
    )

    genotype_column <- ncol(sample_df)

    sample_het <- sample_df[
        sample_df[[genotype_column]] == "0/1",
        ,
        drop = FALSE
    ]

    colnames(sample_het)[genotype_column] <- "sample"

    het_list[[i]] <- sample_het
}


# ------------------------------------------------------------------------------
# 5. Combine individuals
# ------------------------------------------------------------------------------

het_get <- bind_rows(het_list)

message(
    "Total heterozygous variant observations: ",
    nrow(het_get)
)

write.table(
    het_get,
    file = file.path(
        CONFIG$output_dir,
        "TOPMed_heterozygous_PTC_SNVs.tsv"
    ),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
)


# ------------------------------------------------------------------------------
# 6. Construct standardized variant fields
# ------------------------------------------------------------------------------

topmed <- het_get

topmed$CHROM <- topmed$contig
topmed$POS   <- topmed$position

topmed$REF_ALLELE <- vapply(
    strsplit(topmed$variantID, "_"),
    `[`,
    character(1),
    3
)

topmed$ALT_ALLELE <- vapply(
    strsplit(topmed$variantID, "_"),
    `[`,
    character(1),
    4
)


# ------------------------------------------------------------------------------
# 7. Retain SNVs only
# ------------------------------------------------------------------------------

topmed <- topmed[
    nchar(topmed$REF_ALLELE) == 1 &
    nchar(topmed$ALT_ALLELE) == 1,
]

message(
    "Heterozygous SNV observations retained: ",
    nrow(topmed)
)


# ------------------------------------------------------------------------------
# 8. Standardized variant identifiers
# ------------------------------------------------------------------------------

topmed$CHROM <- ifelse(
    grepl("^chr", topmed$CHROM),
    topmed$CHROM,
    paste0("chr", topmed$CHROM)
)

topmed$key <- with(
    topmed,
    paste0(
        CHROM,
        ":",
        POS,
        "_",
        REF_ALLELE,
        ">",
        ALT_ALLELE
    )
)

topmed$variantID_standardized <- with(
    topmed,
    paste(
        CHROM,
        POS,
        REF_ALLELE,
        ALT_ALLELE,
        sep = "_"
    )
)


# ------------------------------------------------------------------------------
# 9. Prepare ANNOVAR VCF-like input
# ------------------------------------------------------------------------------

annovar_vcf <- data.frame(
    CHROM  = topmed$CHROM,
    START  = topmed$POS,
    END    = topmed$POS,
    REF    = topmed$REF_ALLELE,
    ALT    = topmed$ALT_ALLELE,
    ID     = ".",
    ZYG    = "Het",
    FILTER = "PASS",
    stringsAsFactors = FALSE
)

annovar_vcf_file <- file.path(
    CONFIG$output_dir,
    "TOPMed_PTC_SNV.vcf"
)

write.table(
    annovar_vcf,
    annovar_vcf_file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
)


# ------------------------------------------------------------------------------
# 10. Convert to ANNOVAR input
# ------------------------------------------------------------------------------

convert2annovar <- file.path(
    CONFIG$annovar_dir,
    "convert2annovar.pl"
)

annovar_avinput <- file.path(
    CONFIG$output_dir,
    "TOPMed_PTC_SNV.avinput"
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
# 11. ANNOVAR ensGene annotation
# ------------------------------------------------------------------------------

annotate_variation <- file.path(
    CONFIG$annovar_dir,
    "annotate_variation.pl"
)

annovar_output <- file.path(
    CONFIG$output_dir,
    "TOPMed_PTC_SNV_ensGene"
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
    shQuote(annovar_output),
    shQuote(annovar_avinput),
    shQuote(CONFIG$annovar_db)
)

status <- system(annotation_command)

if (status != 0) {
    stop("ANNOVAR gene annotation failed.")
}


message("TOPMed extraction and ANNOVAR preparation complete.")
