############################################################
# 02_LMNA_add_conservation_and_predict.R
# Merge conservation data and export exactly 8 TrunKitten features.
# This script prepares input; prediction is run separately in Python.
############################################################
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript 02_LMNA_add_conservation_and_predict.R <LMNA_RDS> <CONSERVATION_CSV> <OUTPUT_TSV>")
}
LMNA_FILE <- args[1]
CONSERVATION_FILE <- args[2]
OUTPUT_FILE <- args[3]

lmna_snv <- readRDS(LMNA_FILE)
conservation <- read.csv(CONSERVATION_FILE, check.names = FALSE)

required_conservation <- c("variant_id", "phastcons_new3utr_first200_median")
if (!all(required_conservation %in% names(conservation))) {
  stop("Conservation input missing: ", paste(setdiff(required_conservation, names(conservation)), collapse = ", "))
}
if (anyDuplicated(conservation$variant_id)) {
  stop("Conservation input has duplicate variant_id entries; resolve before joining.")
}

lmna_snv_final <- lmna_snv %>%
  left_join(conservation %>% select(all_of(required_conservation)), by = "variant_id")

features <- c(
  "last.EJC", "relativePTClocation", "half_life_PC1",
  "cdsseqs_AU_content", "mut.exon", "cdsseqs_UC_content",
  "phastcons_new3utr_first200_median", "MedianExpression_log2"
)
required <- c("variant_id", features)
missing <- setdiff(required, names(lmna_snv_final))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
if (anyDuplicated(lmna_snv_final$variant_id)) stop("Duplicate variant IDs in LMNA input")
if (anyNA(lmna_snv_final[, required])) {
  stop("Missing values in required model input; inspect annotation before prediction")
}

trunkitten_input_with_id <- lmna_snv_final %>% select(all_of(required))
dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
write.table(trunkitten_input_with_id, OUTPUT_FILE, sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote ", nrow(trunkitten_input_with_id), " variants and ", length(features), " features to ", OUTPUT_FILE)
