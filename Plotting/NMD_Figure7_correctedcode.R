################################################################################
# NMD escape prediction — Panels A-G
#   gnomAD vs ClinVar P/LP vs ClinVar VUS vs GREGoR
#
# Inputs (place in working directory or update paths below):
#   gnomAD_predictions_corrected.csv
#   clinvar_predictions_corrected_PLP_Oct2026.xlsx  (sheet "in", uses clinsig_group: "PLP" / "VUS")
#   gregor_predictions_corrected.csv
#   cds_length_by_gene.csv           (CDS length per gene, for panels D-G)
#
# Required columns: escape_probability, predicted_class (0/1),
#                   predicted_label ("NMD"/"escape"), CHROM, hgnc_symbol,
#                   threshold_used
################################################################################

library(tidyverse)
library(patchwork)
library(scales)
library(readxl)

# ── Palette ───────────────────────────────────────────────────────────────────
COL <- c("gnomAD"       = "#1E9442",   # dark green
         "ClinVar P/LP" = "#7B2D8E",   # dark purple
         "ClinVar VUS"  = "#E879B5",   # pink
         "GREGoR"       = "#E67E22")   # orange
COL_LIGHT <- c("gnomAD"       = "#8ED8A4",
               "ClinVar P/LP" = "#C9A5D9",
               "ClinVar VUS"  = "#F5C5DD",
               "GREGoR"       = "#F5C49A")
DATASETS <- c("gnomAD", "ClinVar P/LP", "ClinVar VUS", "GREGoR")

# ── Read data ─────────────────────────────────────────────────────────────────
gnom <- read_csv("gnomAD_predictions_corrected.csv", show_col_types = FALSE) %>%
  mutate(dataset = "gnomAD")

clin_full <- read_excel("clinvar_predictions_corrected_PLP_Oct2026.xlsx", sheet = "in",
                        .name_repair = "unique_quiet")
plp <- clin_full %>% filter(clinsig_group == "PLP") %>% mutate(dataset = "ClinVar P/LP")
vus <- clin_full %>% filter(clinsig_group == "VUS") %>% mutate(dataset = "ClinVar VUS")

greg <- read_csv("gregor_predictions_corrected.csv", show_col_types = FALSE) %>%
  mutate(dataset = "GREGoR")

cols_keep <- c("dataset", "hgnc_symbol", "CHROM", "escape_probability",
               "predicted_class", "predicted_label", "threshold_used")
both <- bind_rows(
  gnom %>% dplyr::select(all_of(cols_keep)),
  plp  %>% dplyr::select(all_of(cols_keep)),
  vus  %>% dplyr::select(all_of(cols_keep)),
  greg %>% dplyr::select(all_of(cols_keep))
) %>% mutate(dataset = factor(dataset, levels = DATASETS))

# Excel stores 15 significant digits, so the ClinVar threshold differs from the
# CSV value in the 16th decimal; round before checking they agree.
TAU <- unique(round(both$threshold_used, 6))
stopifnot(length(TAU) == 1)

# ── Theme ─────────────────────────────────────────────────────────────────────
FONT <- 17
fig_theme <- theme_minimal(base_size = FONT) + theme(
  plot.title       = element_text(face = "bold", size = FONT + 4,
                                  margin = margin(b = 8)),
  axis.title       = element_text(face = "bold", size = FONT),
  axis.text        = element_text(face = "bold", size = FONT - 2, color = "grey15"),
  axis.line        = element_line(color = "grey15", linewidth = 0.9),
  axis.ticks       = element_line(color = "grey25", linewidth = 0.7),
  axis.ticks.length = unit(0.22, "cm"),
  panel.grid.major = element_line(color = "grey92", linewidth = 0.35),
  panel.grid.minor = element_blank(),
  legend.title     = element_blank(),
  legend.text      = element_text(face = "bold", size = FONT - 4),
  legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA),
  plot.tag         = element_text(face = "bold", size = FONT + 12)
)

# ══════════════════════════════════════════════════════════════════════════════
# PANEL A — overlapping density histograms + tau line
# ══════════════════════════════════════════════════════════════════════════════
pA <- ggplot(both, aes(x = escape_probability, fill = dataset, color = dataset)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 50, alpha = 0.45, position = "identity",
                 linewidth = 0.3) +
  geom_vline(xintercept = TAU, linetype = "dashed",
             color = "grey25", linewidth = 1.0) +
  annotate("text", x = TAU + 0.02, y = Inf, vjust = 4,
           label = sprintf("tau = %.2f", TAU),
           fontface = "bold.italic", size = 4.3, hjust = 0, color = "grey20") +
  scale_fill_manual(values  = COL) +
  scale_color_manual(values = COL) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Escape probability distribution (normalized)",
       x = "Escape probability", y = "Density", tag = "A") +
  fig_theme +
  theme(legend.position = c(0.88, 0.78))

# ══════════════════════════════════════════════════════════════════════════════
# PANEL B — class composition (4 bars per group)
# ══════════════════════════════════════════════════════════════════════════════
panel_b_df <- both %>%
  count(dataset, predicted_label) %>%
  group_by(dataset) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  mutate(class = if_else(predicted_label == "NMD", "NMD-sens.", "Escape"),
         class = factor(class, levels = c("NMD-sens.", "Escape")))

pB <- ggplot(panel_b_df, aes(x = class, y = pct, fill = dataset)) +
  geom_col(position = position_dodge(0.8), width = 0.75,
           color = "grey15", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", pct)),
            position = position_dodge(0.8), vjust = -0.4,
            fontface = "bold", size = 3.4) +
  scale_fill_manual(values = COL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), limits = c(0, 80)) +
  labs(title = "Class composition", x = NULL, y = "Variants (%)", tag = "B") +
  fig_theme +
  theme(legend.position = c(0.88, 0.78))

# ══════════════════════════════════════════════════════════════════════════════
# PANEL C — escape fraction by chromosome
# ══════════════════════════════════════════════════════════════════════════════
chrom_levels <- paste0("chr", 1:22)
panel_c_df <- both %>%
  filter(CHROM %in% chrom_levels) %>%
  mutate(CHROM = factor(CHROM, levels = chrom_levels)) %>%
  group_by(dataset, CHROM) %>%
  summarise(pct_escape = 100 * mean(predicted_label == "escape"),
            .groups = "drop")

mean_lines <- both %>% group_by(dataset) %>%
  summarise(mean_pct = 100 * mean(predicted_label == "escape"), .groups = "drop")

pC <- ggplot(panel_c_df, aes(x = CHROM, y = pct_escape, fill = dataset)) +
  geom_col(position = position_dodge(0.85), width = 0.78,
           color = "grey15", linewidth = 0.3) +
  geom_hline(data = mean_lines, aes(yintercept = mean_pct, color = dataset),
             linetype = "dotted", linewidth = 0.9, show.legend = FALSE) +
  scale_fill_manual(values  = COL) +
  scale_color_manual(values = COL) +
  scale_x_discrete(labels = 1:22) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Escape fraction by chromosome  (dotted lines = group mean)",
       x = "Chromosome", y = "Escape calls (%)", tag = "C") +
  fig_theme +
  theme(legend.position = c(0.06, 0.88),
        legend.direction = "horizontal")

# ══════════════════════════════════════════════════════════════════════════════
# PANELS D-G — variant density per gene = variants / CDS length (per kb)
#   Stacked NMD-sensitive / Escape. Top 15 genes ranked by density among genes
#   with at least `min_n` variants (stops 1-variant short genes dominating).
# ══════════════════════════════════════════════════════════════════════════════

# ── CDS length per gene ──────────────────────────────────────────────────────
#    cds_length_by_gene.csv was built from gnomAD v4.1 constraint metrics
#    (MANE Select CDS; Ensembl canonical / longest CDS as fallback).
#    Keep it in the working directory. The biomaRt block below only runs if
#    the file is missing (requires BiocManager::install("biomaRt")).
CDS_FILE <- "cds_length_by_gene.csv"

if (!file.exists(CDS_FILE)) {
  library(biomaRt)
  genes <- unique(na.omit(c(gnom$hgnc_symbol, plp$hgnc_symbol,
                            vus$hgnc_symbol,  greg$hgnc_symbol)))
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

  # 1) transcripts per gene + canonical flag   (feature page)
  tx <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id",
                             "transcript_biotype", "transcript_is_canonical"),
              filters = "hgnc_symbol", values = genes, mart = mart)
  tx <- tx %>% filter(transcript_biotype == "protein_coding")

  # 2) CDS length per transcript              (structure page)
  cds <- getBM(attributes = c("ensembl_transcript_id", "cds_length"),
               filters = "ensembl_transcript_id",
               values = unique(tx$ensembl_transcript_id), mart = mart) %>%
    distinct() %>% filter(!is.na(cds_length))

  # 3) canonical transcript CDS; fall back to longest CDS if no canonical
  cds_len <- tx %>%
    inner_join(cds, by = "ensembl_transcript_id") %>%
    mutate(is_canon = !is.na(transcript_is_canonical) &
                      transcript_is_canonical == 1) %>%
    group_by(hgnc_symbol) %>%
    arrange(desc(is_canon), desc(cds_length), .by_group = TRUE) %>%
    slice(1) %>% ungroup() %>%
    dplyr::select(hgnc_symbol, ensembl_transcript_id, cds_length)

  write_csv(cds_len, CDS_FILE)
}
cds_len <- read_csv(CDS_FILE, show_col_types = FALSE)

make_top_genes_plot <- function(df, dataset_name, col_dark, col_light, tag,
                                min_n = 5) {
  df <- df %>% filter(!is.na(hgnc_symbol), hgnc_symbol != "")

  gene_tot <- df %>%
    count(hgnc_symbol, name = "total") %>%
    inner_join(cds_len %>% dplyr::select(hgnc_symbol, cds_length),
               by = "hgnc_symbol") %>%
    mutate(cds_kb  = cds_length / 1000,
           density = total / cds_kb)

  n_missing <- n_distinct(df$hgnc_symbol) - nrow(gene_tot)
  if (n_missing > 0)
    message(sprintf("%s: %d genes without CDS length were dropped",
                    dataset_name, n_missing))

  top15 <- gene_tot %>%
    filter(total >= min_n) %>%
    slice_max(density, n = 15, with_ties = FALSE)

  plot_df <- df %>%
    filter(hgnc_symbol %in% top15$hgnc_symbol) %>%
    count(hgnc_symbol, predicted_label) %>%
    left_join(top15, by = "hgnc_symbol") %>%
    mutate(seg_density = n / cds_kb,
           class = factor(if_else(predicted_label == "NMD",
                                  "NMD-sensitive", "Escape"),
                          levels = c("NMD-sensitive", "Escape")),
           hgnc_symbol = fct_reorder(hgnc_symbol, density))

  labels_df <- top15 %>%
    mutate(hgnc_symbol = factor(hgnc_symbol,
                                levels = levels(plot_df$hgnc_symbol)))

  ggplot(plot_df, aes(x = seg_density, y = hgnc_symbol, fill = class)) +
    geom_col(color = "grey15", linewidth = 0.3) +
    geom_text(data = labels_df,
              aes(x = density, y = hgnc_symbol,
                  label = sprintf("%.2f  (n=%d)", density, total)),
              inherit.aes = FALSE, hjust = -0.06,
              fontface = "bold", size = 3.2, color = "grey25") +
    scale_fill_manual(values = c("NMD-sensitive" = unname(col_dark),
                                 "Escape"        = unname(col_light))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
    labs(title = sprintf("Top 15 genes by variant density  (%s)", dataset_name),
         subtitle = sprintf("genes with \u2265 %d variants", min_n),
         x = "Variants per kb of CDS", y = NULL, tag = tag) +
    fig_theme +
    theme(axis.text.y     = element_text(face = "bold.italic", size = FONT),
          plot.subtitle   = element_text(size = FONT - 4, color = "grey35"),
          legend.position = c(0.85, 0.12),
          legend.background = element_rect(fill = alpha("white", 0.9),
                                           color = "grey80", linewidth = 0.3))
}

pD <- make_top_genes_plot(gnom, "gnomAD",       COL[["gnomAD"]],
                          COL_LIGHT[["gnomAD"]],       "D")
pE <- make_top_genes_plot(plp,  "ClinVar P/LP", COL[["ClinVar P/LP"]],
                          COL_LIGHT[["ClinVar P/LP"]], "E")
pF <- make_top_genes_plot(vus,  "ClinVar VUS",  COL[["ClinVar VUS"]],
                          COL_LIGHT[["ClinVar VUS"]],  "F")
pG <- make_top_genes_plot(greg, "GREGoR",       COL[["GREGoR"]],
                          COL_LIGHT[["GREGoR"]],       "G",
                          min_n = 3)   # GREGoR is small: only 17 genes have >= 5 variants

# ══════════════════════════════════════════════════════════════════════════════
# COMBINE
#   Row 1: A | B
#   Row 2: C  (full width)
#   Row 3: D | E
#   Row 4: F | G
# ══════════════════════════════════════════════════════════════════════════════
combined <- (pA | pB) /
             pC        /
            (pD | pE) /
            (pF | pG) +
  plot_layout(heights = c(1.0, 0.9, 1.4, 1.4))

ggsave("NMD_escape_4datasets_panels_A_to_G.pdf", combined,
       width = 17, height = 22, dpi = 300)
ggsave("NMD_escape_4datasets_panels_A_to_G.png", combined,
       width = 17, height = 22, dpi = 300)

cat("Saved: NMD_escape_4datasets_panels_A_to_G.pdf / .png\n")
print(combined)
