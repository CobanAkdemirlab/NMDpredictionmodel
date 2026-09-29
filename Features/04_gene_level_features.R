# ==============================================================================
# Gene-level feature annotation
#
# Purpose:
#   Add gene-level constraint, expression, and mRNA stability features
#   to the annotated PTC variant datasets.
#
# Shared across:
#   TOPMed, gnomAD, ClinVar, GREGoR
# ==============================================================================

library(dplyr)
library(readxl)

CONFIG <- list(
    # gnomAD v4.1 gene constraint metrics
    #
    # Source:
    # https://gnomad.broadinstitute.org/downloads
    #
    # Download:
    # gnomAD v4.1 -> Constraint -> Constraint metrics
    #
    # File:
    # gnomad.v4.1.constraint_metrics.tsv
    #
    # The full gnomAD file is not distributed with this repository.
    # Download it from gnomAD and update the path below.
    lof_metrics_file =
        "/path/to/gnomad.v4.1.constraint_metrics.tsv",

    # GTEx v8 median gene expression by tissue
    #
    # Source:
    # https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
    #
    # Download:
    # GTEx Analysis V8 -> RNA-Seq ->
    # Median gene-level TPM by tissue
    #
    # File:
    # GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
    #
    # GTEx reports that these tissue-level median TPM values were calculated
    # from GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz.
    
    gtex_expression_file =
        "/path/to/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.txt",

    # Agarwal & Kelley (2022) supplementary human mRNA half-life data
    #
    # Source:
    # Agarwal V, Kelley DR. The genetic and biochemical determinants of
    # mRNA degradation rates in mammals.
    # Genome Biol. 2022;23:245.
    # doi: 10.1186/s13059-022-02811-x
    #
    # File:
    # 13059_2022_2811_MOESM3_ESM.xlsx
    half_life_file =
        "/path/to/13059_2022_2811_MOESM3_ESM.xlsx",
    
    # Transcript-to-gene mapping used to assign stable Ensembl gene IDs
    # to the annotated variant transcripts.
    canonical_gene_map =
        "/path/to/canonical_transcript_gene_map.tsv"
)

# ------------------------------------------------------------------------------
# 1. Load annotated variants
# ------------------------------------------------------------------------------

variants <- readRDS(
    "/path/to/PTC_annotated_variants.rds"
)


# ------------------------------------------------------------------------------
# 2. Add stable Ensembl gene mapping
# ------------------------------------------------------------------------------

gene_map <- read.delim(
    CONFIG$canonical_gene_map,
    stringsAsFactors = FALSE
)

variants <- variants %>%
    left_join(
        gene_map %>%
            select(
                ensembl_transcript_id,
                ensembl_gene_id,
                hgnc_symbol
            ),
        by = c(
            "txnames" = "ensembl_transcript_id"
        )
    )


# ------------------------------------------------------------------------------
# 3. gnomAD gene constraint
# ------------------------------------------------------------------------------

lof_metrics <- read.table(
    CONFIG$lof_metrics_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

# Retain one canonical transcript per gene
lof_sub <- lof_metrics %>%
    filter(
        canonical == 'true'
    ) %>%
    select(
        gene_id,
        pLI          = lof.pLI,
        oe_lof_upper = lof.oe_ci.upper
    )

# QC: canonical gene IDs should be unique
stopifnot(
    anyDuplicated(lof_sub$gene_id) == 0
)

# Merge gene-level constraint metrics
variants <- variants %>%
    left_join(
        lof_sub,
        by = c(
            "ensembl_gene_id" = "gene_id"
        ),
        relationship = "many-to-one"
    )

# Categorization 
# pli

ind <- which(df$pLI<0.35 )
df$pLI.cat <- rep('NA',nrow(df))
df$pLI.cat[ind] <- 'highly tolerant (pLI<0.35)'


ind <- which(df$pLI>=0.35 & df$pLI<0.65 )
df$pLI.cat[ind] <- 'medium tolerant (0.35=<pLI<0.65)'


ind <- which(df$pLI>=0.65 )
df$pLI.cat[ind] <- 'highly intolerant (pLI>=0.65)'

#level of intolerance
new_topmed$loeuf_cat = cut(new_topmed$oe_lof_upper, breaks = c(0,0.2,0.6,2), 
                           labels = c('Very LoF-intolerant (0-0.2)','Moderately_constrained (0.2-0.6)','LoF_tolerant (0.6-2)')) #level of intolerance

# ------------------------------------------------------------------------------
# 4. GTEx expression
# ------------------------------------------------------------------------------

expression <- read.table(
    CONFIG$gtex_expression_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

expression$gene_id <- sub(
    "\\..*",
    "",
    expression$Name
)

tissue_cols <- setdiff(
    colnames(expression),
    c(
        "Name",
        "Description",
        "gene_id"
    )
)

expression$MedianExpression <- apply(
    expression[, tissue_cols],
    1,
    median,
    na.rm = TRUE
)

expression_sub <- expression %>%
    select(
        gene_id,
        MedianExpression,
        Whole.Blood
    )

variants <- variants %>%
    left_join(
        expression_sub,
        by = c(
            "ensembl_gene_id" = "gene_id"
        )
    ) %>%
    mutate(
        MedianExpression_log2 =
            log2(MedianExpression + 1)
    )


# ------------------------------------------------------------------------------
# 5. mRNA half-life PC1
# ------------------------------------------------------------------------------
# Source:
# Agarwal V, Kelley DR. The genetic and biochemical determinants of
# mRNA degradation rates in mammals. Genome Biology. 2022;23:245.
# https://doi.org/10.1186/s13059-022-02811-x
#
# The study compiled transcriptome-wide mRNA decay measurements across
# multiple human datasets and derived a consensus, cell-type-agnostic
# measure of mRNA half-life using the first principal component (PC1).
#
# The human consensus half-life values were obtained from the study's
# supplementary data:
#
#   13059_2022_2811_MOESM3_ESM.xlsx
#
# Feature used in this study:
#   half_life_PC1

half_life <- read_excel(
    CONFIG$half_life_file,
    sheet = "human",
    skip = 1
)

half_life_sub <- half_life %>%
    select(
        ensembl_gene_id =
            `Ensembl Gene Id`,
        half_life_PC1 =
            `half-life (PC1)`
    )

variants <- variants %>%
    left_join(
        half_life_sub,
        by = "ensembl_gene_id"
    )

