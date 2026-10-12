############################################################
# 03_LMNA_SGE_TrunKitten_validation.R
# External validation against Cortazar LMNA SGE measurements.
# Spearman correlation: (1 - predicted escape probability) vs
# region-normalized experimental NMD measurement.
############################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript 03_LMNA_SGE_TrunKitten_validation.R <PTC_FILE> <SNV_FILE> <LMNA_RDS> <PREDICTIONS_TSV> <OUTPUT_DIR>")
}
PTC_FILE <- args[1]
SNV_FILE <- args[2]
LMNA_FILE <- args[3]
PREDICTION_FILE <- args[4]
OUTPUT_DIR <- args[5]
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

PTC <- read.table(PTC_FILE, header = TRUE)
SNV <- read.table(SNV_FILE, header = FALSE)
if (ncol(SNV) != 3) stop("Expected three columns in SNV experimental input")
names(SNV) <- c("codon", "mutation", "NMD")
lmna_snv <- readRDS(LMNA_FILE)
pred <- read.delim(PREDICTION_FILE, check.names = FALSE)

get_region_median <- function(data, min_codon, max_codon) {
  data %>%
    mutate(codon = as.numeric(codon), NMD = as.numeric(NMD)) %>%
    filter(codon >= min_codon, codon <= max_codon) %>%
    group_by(codon, mutation) %>%
    summarise(SNVexp_NMD = mean(NMD), .groups = "drop") %>%
    summarise(SNVmedian_NMD = median(SNVexp_NMD)) %>%
    pull(SNVmedian_NMD)
}

# Preserve the regional normalization and codon assignments from the
# original Cortazar LMNA validation workflow.
regions <- data.frame(
  snv_start = c(1,28,59,89,214,247,434,537,567),
  snv_end   = c(27,58,88,118,246,270,460,566,656),
  ptc_start = c(1,24,59,87,214,240,434,537,567),
  ptc_end   = c(23,58,86,118,239,270,460,566,656)
)
region_tables <- lapply(seq_len(nrow(regions)), function(i) {
  r <- regions[i, ]
  med <- get_region_median(SNV, r$snv_start, r$snv_end)
  tibble(codon = seq.int(r$ptc_start, r$ptc_end), SNVmedian_NMD = med)
})
dfx_SNVs_medians <- bind_rows(region_tables)

PTC_norm <- PTC %>%
  select(codon, mutation, mean_val, std_val) %>%
  mutate(codon = as.numeric(codon), mean_val = as.numeric(mean_val), std_val = as.numeric(std_val)) %>%
  filter(!is.na(std_val), std_val <= 0.2) %>%
  left_join(dfx_SNVs_medians, by = "codon") %>%
  mutate(norm_mean = mean_val - SNVmedian_NMD)

if (anyDuplicated(PTC_norm[c("codon", "mutation")])) {
  stop("Duplicate PTC codon/mutation keys; joining would duplicate validation variants")
}
if (anyDuplicated(lmna_snv$variant_id)) stop("Duplicate LMNA variant_id entries")
if (!all(c("variant_id", "escape_prob", "escape_pred_at_youden", "threshold_used") %in% names(pred))) {
  stop("Prediction TSV lacks expected columns; inspect predict.py output")
}
if (anyDuplicated(pred$variant_id)) stop("Duplicate prediction variant_id entries")

lmna_validation <- lmna_snv %>%
  left_join(PTC_norm %>% select(codon, stop_codon = mutation, SNVmedian_NMD, norm_mean),
            by = c("codon", "stop_codon")) %>%
  left_join(pred %>% select(variant_id, escape_prob, escape_pred_at_youden, threshold_used),
            by = "variant_id")
if (anyNA(lmna_validation$escape_prob)) stop("Some LMNA variants lack model predictions")
if (any(lmna_validation$escape_prob < 0 | lmna_validation$escape_prob > 1)) {
  stop("escape_prob values must be within [0,1]")
}
message("LMNA variants: ", nrow(lmna_validation))
message("Variants lacking normalized experimental NMD: ", sum(is.na(lmna_validation$norm_mean)))
message("Threshold(s) in predictions: ", paste(unique(lmna_validation$threshold_used), collapse = ", "))

lmna_validation_final <- lmna_validation %>%
  filter(!is.na(norm_mean)) %>%
  mutate(predicted_NMD = 1 - escape_prob)
if (nrow(lmna_validation_final) < 3) stop("Insufficient paired observations")

result <- cor.test(lmna_validation_final$predicted_NMD,
                   lmna_validation_final$norm_mean,
                   method = "spearman", exact = FALSE)
print(result)
summary_df <- data.frame(n = nrow(lmna_validation_final),
                         spearman_rho = unname(result$estimate), p_value = result$p.value)
write.csv(summary_df, file.path(OUTPUT_DIR, "LMNA_SGE_validation_statistics.csv"), row.names = FALSE)
write.csv(lmna_validation_final, file.path(OUTPUT_DIR, "LMNA_SGE_validation_final.csv"), row.names = FALSE)

p <- ggplot(lmna_validation_final,
            aes(x = predicted_NMD, y = norm_mean, fill = factor(mut.exon))) +
  geom_point(shape = 21, color = "black", size = 2.7, stroke = 0.5) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
  theme_classic(base_size = 14) +
  labs(x = "TrunKitten-predicted NMD probability",
       y = "Region-normalized experimental NMD measurement",
       fill = "LMNA exon")
ggsave(file.path(OUTPUT_DIR, "LMNA_SGE_TrunKitten_validation.png"),
       plot = p, width = 7, height = 5, dpi = 300)
message("Validation outputs written to: ", OUTPUT_DIR)
