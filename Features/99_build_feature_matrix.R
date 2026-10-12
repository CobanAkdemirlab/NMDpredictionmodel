# ==============================================================================
# 99_build_feature_matrix.R
# Assemble core feature modules 01-07 for TOPMed, gnomAD, ClinVar, GREGoR.
# Modules 08-13 are merged downstream in the TrunCat notebook.
# Run from the repository root, e.g.:
#   NMD_DATASET=ClinVar Rscript Features/99_build_feature_matrix.R
# Optional: override individual paths using the environment variables below.
# ==============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

DATASET <- Sys.getenv("NMD_DATASET", "TOPMed")
stopifnot(DATASET %in% c("TOPMed", "gnomAD", "ClinVar", "GREGoR"))

output_dir <- file.path("output", DATASET)
feature_dir <- file.path(output_dir, "features")

# Environment-variable overrides allow the exact output paths used by each
# feature script without editing this shared assembler.
path_setting <- function(env, default) {
    value <- Sys.getenv(env, unset = "")
    if (nzchar(value)) value else default
}

config <- list(
    base_variants = path_setting("NMD_BASE_VARIANTS", file.path(output_dir, paste0(DATASET, "_filtered_variants.rds"))),
    gencode_features = path_setting("NMD_01A_RDS", file.path(feature_dir, "gencode_v26_sequence_features.rds")),
    intron_features = Sys.getenv("NMD_01B_CSV", unset = ""),
    ptbp1_features = path_setting("NMD_02_RDS", file.path(feature_dir, "02_PTBP1_binding_features.rds")),
    amino_acid_features = path_setting("NMD_03_RDS", file.path(feature_dir, "PTC_amino_acid_features.rds")),
    gene_features = path_setting("NMD_04_RDS", file.path(feature_dir, "04_gene_level_features.rds")),
    ptc_geometry_features = path_setting("NMD_05_RDS", file.path(feature_dir, "PTC_geometry_features.rds")),
    variant_scores = path_setting("NMD_06_RDS", file.path(feature_dir, "variant_annotation_scores.rds")),
    motif_features = path_setting("NMD_07_CSV", file.path(feature_dir, "07_FIMO_motif_features.csv")),
    output = path_setting("NMD_CORE_OUTPUT", file.path(output_dir, paste0(DATASET, "_core_features.csv")) )
)

required_paths <- config[c("base_variants", "gencode_features", "ptbp1_features",
                           "amino_acid_features", "gene_features", "ptc_geometry_features",
                           "variant_scores", "motif_features")]
missing <- names(required_paths)[!file.exists(unlist(required_paths))]
if (length(missing)) {
    stop("Missing required inputs (set the corresponding NMD_* path or generate the output):\n",
         paste(sprintf("  %s: %s", missing, unlist(required_paths[missing])), collapse = "\n"),
         call. = FALSE)
}
if (nzchar(config$intron_features) && !file.exists(config$intron_features)) {
    stop("NMD_01B_CSV not found: ", config$intron_features)
}

read_rds_table <- function(path, label) {
    x <- readRDS(path)
    if (!is.data.frame(x)) stop(label, " must be a data.frame")
    as_tibble(x)
}

validate_key <- function(x, label, unique_required = TRUE) {
    if (!"key" %in% names(x)) stop(label, " is missing the `key` column")
    if (anyNA(x$key) || any(!nzchar(as.character(x$key)))) {
        stop(label, " contains missing or blank keys")
    }
    if (unique_required && anyDuplicated(x$key)) {
        stop(label, " contains duplicate keys (feature tables must have one row per genomic key)")
    }
    message(label, ": ", nrow(x), " rows; ", ncol(x), " columns")
    invisible(x)
}

base <- read_rds_table(config$base_variants, "Base variants")
# The base may legitimately contain multiple ClinVar records for one genomic
# variant. A many-to-one join preserves those records without multiplying rows.
validate_key(base, "Base variants", unique_required = FALSE)

features <- list(
    GENCODE = read_rds_table(config$gencode_features, "GENCODE"),
    PTBP1 = read_rds_table(config$ptbp1_features, "PTBP1"),
    amino_acid = read_rds_table(config$amino_acid_features, "Amino acid"),
    gene = read_rds_table(config$gene_features, "Gene"),
    PTC_geometry = read_rds_table(config$ptc_geometry_features, "PTC geometry"),
    variant_scores = read_rds_table(config$variant_scores, "Variant scores"),
    motifs = read_csv(config$motif_features, show_col_types = FALSE)
)

# Module 04 currently annotates an entire variant table. Keep its gene-level
# outputs only, so unrelated columns cannot overwrite core annotations.
gene_cols <- c("key", "ensembl_gene_id", "hgnc_symbol", "pLI", "oe_lof_upper",
               "pLI.cat", "loeuf_cat", "MedianExpression", "Whole.Blood",
               "MedianExpression_log2", "half_life_PC1")
features$gene <- features$gene %>% select(any_of(gene_cols))

# 01b edits three intron columns in a CSV rather than writing a separate RDS.
# Supply that CSV with NMD_01B_CSV. If the base already has those columns,
# this step replaces their values by key; otherwise it adds them.
intron_cols <- c("cdsseq.introns", "threeUTR.introns", "fiveUTR.introns")
if (nzchar(config$intron_features)) {
    intron_table <- read_csv(config$intron_features, show_col_types = FALSE)
    if (!"key" %in% names(intron_table)) {
        stop("01b CSV requires a `key` column for a safe merge")
    }
    if (!all(intron_cols %in% names(intron_table))) {
        stop("01b CSV is missing required intron columns: ",
             paste(setdiff(intron_cols, names(intron_table)), collapse = ", "))
    }
    features <- append(features, list(introns = intron_table %>% select(key, all_of(intron_cols))), after = 1)
} else if (!all(intron_cols %in% names(base))) {
    stop("Intron columns absent from base. Supply 01b's completed CSV via NMD_01B_CSV")
}

for (label in names(features)) validate_key(features[[label]], label)

# For most modules, existing metadata columns are retained from the base.
# For 01b and 04, explicitly replace existing columns with the updated values.
feature_matrix <- base
for (label in names(features)) {
    x <- features[[label]]
    overlapping <- setdiff(intersect(names(feature_matrix), names(x)), "key")
    if (label %in% c("introns", "gene")) {
        feature_matrix <- feature_matrix %>% select(-any_of(overlapping))
    } else if (length(overlapping)) {
        message(label, ": ignoring already-present columns: ", paste(overlapping, collapse = ", "))
        x <- x %>% select(-all_of(overlapping))
    }
    before <- nrow(feature_matrix)
    feature_matrix <- feature_matrix %>% left_join(x, by = "key", relationship = "many-to-one")
    if (nrow(feature_matrix) != before) stop("Row count changed while joining ", label)
    message("Merged ", label, ": ", nrow(feature_matrix), " rows x ", ncol(feature_matrix), " columns")
}

# Retain all base rows and their original order; duplicate genomic keys are
# allowed only in the base (e.g., two ClinVar accessions at one locus).
if (nrow(feature_matrix) != nrow(base)) stop("Base row count not preserved")

missingness <- tibble(
    feature = names(feature_matrix),
    n_missing = vapply(feature_matrix, function(x) sum(is.na(x)), integer(1)),
    pct_missing = vapply(feature_matrix, function(x) round(mean(is.na(x)) * 100, 2), numeric(1))
) %>% arrange(desc(pct_missing))

dir.create(dirname(config$output), recursive = TRUE, showWarnings = FALSE)
write_csv(feature_matrix, config$output, na = "")
qc_path <- file.path(dirname(config$output), paste0(DATASET, "_core_feature_missingness.csv"))
write_csv(missingness, qc_path)
message("Saved core table: ", config$output)
message("Saved missingness QC: ", qc_path)
message("Records: ", nrow(feature_matrix), "; unique genomic keys: ", n_distinct(feature_matrix$key))

if (DATASET == "TOPMed") {
    truncat_output <- file.path("Model", "TrunCat", "data", "TOPMed_stopgain.csv")
    dir.create(dirname(truncat_output), recursive = TRUE, showWarnings = FALSE)
    write_csv(feature_matrix, truncat_output, na = "")
    message("Saved TOPMed TrunCat handoff: ", truncat_output)
}
message("Core assembly complete. Merge modules 08-13 downstream.")
