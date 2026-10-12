# ==============================================================================
# 99_build_feature_matrix.R
# Assemble core features 01-07. Specialized features 08+ are merged downstream.
# Run from the repository root:
#   NMD_DATASET=ClinVar Rscript Features/99_build_feature_matrix.R
# Input paths may be overridden with NMD_BASE_RDS, NMD_01_RDS,
# NMD_02_RDS, NMD_03_RDS, NMD_04_RDS, NMD_05_RDS, NMD_06_RDS,
# NMD_07_CSV. No upstream feature script is modified.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

DATASET <- Sys.getenv("NMD_DATASET", "TOPMed")
stopifnot(DATASET %in% c("TOPMed", "gnomAD", "ClinVar", "GREGoR"))

path_default <- function(env, default) {
  value <- Sys.getenv(env, "")
  if (nzchar(value)) value else default
}

if (DATASET == "ClinVar") {
  defaults <- list(
    base = "Annotation/output/ClinVar_PTC_filtered.rds",
    sequence = "Features/output/ClinVar_gencode_v26_sequence_features.rds",
    ptbp1 = "Features/output/ClinVar_PTBP1_3UTR_features.rds",
    amino = "Features/output/ClinVar_PTC_amino_acid_features.rds",
    gene = "Features/output/ClinVar_gene_level_features.rds",
    geometry = "Features/output/ClinVar_PTC_geometry_features.rds",
    scores = "Features/output/ClinVar_variant_annotation_scores.rds",
    motifs = "Features/output/ClinVar_FIMO_results/motif_presence_matrix.csv"
  )
} else {
  # Defaults for other cohorts are conventional locations; override to match
  # the actual files produced by their feature runs.
  d <- file.path("output", DATASET)
  f <- file.path(d, "features")
  defaults <- list(
    base = file.path(d, paste0(DATASET, "_filtered_variants.rds")),
    sequence = file.path(f, "gencode_v26_sequence_features.rds"),
    ptbp1 = file.path(f, "02_PTBP1_binding_features.rds"),
    amino = file.path(f, "PTC_amino_acid_features.rds"),
    gene = file.path(f, "04_gene_level_features.rds"),
    geometry = file.path(f, "PTC_geometry_features.rds"),
    scores = file.path(f, "variant_annotation_scores.rds"),
    motifs = file.path(f, "motif_presence_matrix.csv")
  )
}

config <- list(
  base = path_default("NMD_BASE_RDS", defaults$base),
  sequence = path_default("NMD_01_RDS", defaults$sequence),
  ptbp1 = path_default("NMD_02_RDS", defaults$ptbp1),
  amino = path_default("NMD_03_RDS", defaults$amino),
  gene = path_default("NMD_04_RDS", defaults$gene),
  geometry = path_default("NMD_05_RDS", defaults$geometry),
  scores = path_default("NMD_06_RDS", defaults$scores),
  motifs = path_default("NMD_07_CSV", defaults$motifs),
  output = path_default("NMD_CORE_OUTPUT", file.path("output", DATASET,
                              paste0(DATASET, "_core_features.csv")))
)

for (nm in setdiff(names(config), "output")) {
  if (!file.exists(config[[nm]])) {
    stop("Missing ", nm, " input: ", config[[nm]], call. = FALSE)
  }
}

base <- as.data.frame(readRDS(config$base))
if (!"key" %in% names(base) || anyNA(base$key) || any(!nzchar(as.character(base$key)))) {
  stop("Base table needs a complete `key` column.")
}
if (!"txnames" %in% names(base) || anyNA(base$txnames)) {
  stop("Base table needs a complete `txnames` column for transcript features.")
}
base$.nmd_record_order <- seq_len(nrow(base))
message("Dataset: ", DATASET, "; base records: ", nrow(base),
        "; unique keys: ", dplyr::n_distinct(base$key))

# Collapse duplicate identifiers ONLY when all non-ID values agree.
# This preserves repeated ClinVar records in the base table while preventing
# one-to-many feature joins. Conflicts are reported rather than silently lost.
unique_feature_rows <- function(x, id, label) {
  x <- as.data.frame(x)
  if (!id %in% names(x)) stop(label, " lacks `", id, "`.")
  if (anyNA(x[[id]]) || any(!nzchar(as.character(x[[id]])))) {
    stop(label, " has missing identifiers.")
  }
  if (anyDuplicated(names(x))) stop(label, " has duplicated column names.")
  cols <- setdiff(names(x), id)
  if (anyDuplicated(x[[id]])) {
    groups <- split(seq_len(nrow(x)), as.character(x[[id]]))
    for (idx in groups) {
      if (length(idx) < 2L) next
      for (col in cols) {
        vals <- x[[col]][idx]
        if (length(unique(as.character(vals))) > 1L &&
            !(all(is.na(vals)))) {
          stop(label, " has conflicting `", col, "` values for ",
               id, "=", as.character(x[[id]][idx[1]]),
               ". A key-only join would be ambiguous.", call. = FALSE)
        }
      }
    }
    x <- x[!duplicated(x[[id]]), , drop = FALSE]
  }
  x
}

# Apply one feature family at a time. New columns are added; existing columns
# are filled only where missing. Conflicting nonmissing values are flagged.
merge_family <- function(target, source, id, label) {
  source <- unique_feature_rows(source, id, label)
  # Metadata already in the base table should not be merged as new features.
  # If shared feature columns exist, fill missing values and check conflicts.
  shared <- intersect(setdiff(names(source), id), names(target))
  source_new <- source[, c(id, setdiff(names(source), names(target))), drop = FALSE]
  before <- nrow(target)
  target <- dplyr::left_join(target, source_new, by = id, relationship = "many-to-one")
  if (nrow(target) != before) stop("Row count changed merging ", label)
  if (length(shared)) {
    lookup <- match(as.character(target[[id]]), as.character(source[[id]]))
    for (nm in shared) {
      # Base identifiers/annotation metadata are authoritative; don't overwrite.
      # Introns are features and can be filled when missing.
      if (!nm %in% c("cdsseq.introns", "threeUTR.introns", "fiveUTR.introns")) next
      incoming <- source[[nm]][lookup]
      missing <- is.na(target[[nm]]) | (is.character(target[[nm]]) & target[[nm]] == "")
      target[[nm]][missing] <- incoming[missing]
    }
  }
  message(sprintf("  %-14s join=%-7s feature rows=%d; matrix=%d x %d",
                  label, id, nrow(source), nrow(target), ncol(target)))
  target
}

# 01a already includes the 01b intron columns in the saved sequence table.
# 02 is also transcript-level. All other inputs are variant-keyed.
feature_specs <- list(
  list(name="GENCODE", path=config$sequence, id="txnames", type="rds"),
  list(name="PTBP1", path=config$ptbp1, id="txnames", type="rds"),
  list(name="amino_acid", path=config$amino, id="key", type="rds"),
  list(name="gene", path=config$gene, id="key", type="rds"),
  list(name="PTC_geometry", path=config$geometry, id="key", type="rds"),
  list(name="variant_scores", path=config$scores, id="key", type="rds"),
  list(name="motifs", path=config$motifs, id="key", type="csv")
)

feature_matrix <- base
for (spec in feature_specs) {
  feature <- if (spec$type == "csv") {
    as.data.frame(readr::read_csv(spec$path, show_col_types = FALSE,
                                   name_repair = "minimal", progress = FALSE))
  } else {
    as.data.frame(readRDS(spec$path))
  }
  # FIMO contains fasta_id/txnames metadata. Only key and motif columns are
  # retained; transcript metadata is already present in the base table.
  if (spec$name == "motifs") {
    feature <- feature[, setdiff(names(feature), c("fasta_id", "txnames")), drop=FALSE]
  }
  feature_matrix <- merge_family(feature_matrix, feature, spec$id, spec$name)
}

feature_matrix <- feature_matrix[order(feature_matrix$.nmd_record_order), , drop=FALSE]
feature_matrix$.nmd_record_order <- NULL
if (nrow(feature_matrix) != nrow(base)) stop("Final row count mismatch.")

missingness <- data.frame(
  feature = names(feature_matrix),
  n_missing = vapply(feature_matrix, function(x) sum(is.na(x)), integer(1)),
  pct_missing = vapply(feature_matrix, function(x) round(mean(is.na(x))*100, 2), numeric(1)),
  stringsAsFactors = FALSE
)
missingness <- dplyr::arrange(missingness, dplyr::desc(pct_missing))
dir.create(dirname(config$output), recursive=TRUE, showWarnings=FALSE)
readr::write_csv(feature_matrix, config$output, na="")
qc_path <- file.path(dirname(config$output), paste0(DATASET, "_core_feature_missingness.csv"))
readr::write_csv(missingness, qc_path)
message("Saved: ", config$output)
message("Saved QC: ", qc_path)
message("Final: ", nrow(feature_matrix), " records; ", ncol(feature_matrix),
        " columns; ", dplyr::n_distinct(feature_matrix$key), " unique genomic keys.")

if (DATASET == "TOPMed") {
  handoff <- "Model/TrunCat/data/TOPMed_stopgain.csv"
  dir.create(dirname(handoff), recursive=TRUE, showWarnings=FALSE)
  readr::write_csv(feature_matrix, handoff, na="")
  message("Saved TOPMed TrunCat handoff: ", handoff)
}
message("Core feature assembly complete. Features 08+ remain separate.")
