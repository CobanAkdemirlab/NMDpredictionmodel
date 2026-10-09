################################################################################
#   Supplemetal Figure 12 
#   NMD escape prediction — Panels A-G (capital letters)
#   Restricted to genes that are BOTH:
#     (1) HIGHLY INTOLERANT  (LOEUF < 0.60  AND  pLI >= 0.90)   [v4.1 default]
#     (2) AUTOSOMAL DOMINANT (omim_AD_symbols.csv)
#
# Datasets: gnomAD, ClinVar P/LP, ClinVar VUS, GREGoR
#

#
# Inputs (same directory or update paths below):
#   gnomAD_predictions_corrected.csv
#   clinvar_predictions_corrected_PLP_Oct2026.xlsx
#   gregor_predictions_corrected.csv
#   gnomad_v4_1_constraint_metrics.tsv
#   omim_AD_symbols.csv
#   cds_length_by_gene.csv
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(readxl)
})

# ── User-configurable ────────────────────────────────────────────────────────
GNOM_PATH       <- "gnomAD_predictions_corrected.csv"
CLIN_PATH       <- "clinvar_predictions_corrected_PLP_Oct2026.xlsx"
CLIN_SHEET      <- "in"
GREG_PATH       <- "gregor_predictions_corrected.csv"
CONSTRAINT_PATH <- "gnomad.v4.1.constraint_metrics.tsv"
AD_GENES_PATH   <- "omim_AD_symbols.csv"
CDS_FILE        <- "cds_length_by_gene.csv"

LOEUF_MAX <- 0.60     # v4.1 recommended (use 0.35 for v2.1.1)
PLI_MIN   <- 0.90
MIN_N_DEFAULT <- 5    # min variants per gene for panels D-G (Figure 7: 5; GREGoR 3)

# ── Palette ───────────────────────────────────────────────────────────────────
COL <- c("gnomAD"       = "#1E9442",
         "ClinVar P/LP" = "#7B2D8E",
         "ClinVar VUS"  = "#E879B5",
         "GREGoR"       = "#E67E22")
COL_LIGHT <- c("gnomAD"       = "#8ED8A4",
               "ClinVar P/LP" = "#C9A5D9",
               "ClinVar VUS"  = "#F5C5DD",
               "GREGoR"       = "#F5C49A")
DATASETS <- c("gnomAD", "ClinVar P/LP", "ClinVar VUS", "GREGoR")

# ── Read predictions ─────────────────────────────────────────────────────────
gnom_raw <- read_csv(GNOM_PATH, show_col_types = FALSE)
greg_raw <- read_csv(GREG_PATH, show_col_types = FALSE)
clin_raw <- read_excel(CLIN_PATH, sheet = CLIN_SHEET, .name_repair = "unique_quiet")

stopifnot("clinsig_group" %in% names(clin_raw))
plp_raw <- clin_raw %>% filter(clinsig_group == "PLP")
vus_raw <- clin_raw %>% filter(clinsig_group == "VUS")

cols_keep <- c("hgnc_symbol", "CHROM", "escape_probability",
               "predicted_class", "predicted_label", "threshold_used")

both_all <- bind_rows(
  gnom_raw %>% mutate(group = "gnomAD")       %>% dplyr::select(group, all_of(cols_keep)),
  plp_raw  %>% mutate(group = "ClinVar P/LP") %>% dplyr::select(group, all_of(cols_keep)),
  vus_raw  %>% mutate(group = "ClinVar VUS")  %>% dplyr::select(group, all_of(cols_keep)),
  greg_raw %>% mutate(group = "GREGoR")       %>% dplyr::select(group, all_of(cols_keep))
) %>% mutate(group = factor(group, levels = DATASETS))

# Excel stores 15 significant digits -> round before checking agreement
TAU <- unique(round(both_all$threshold_used, 6))
stopifnot(length(TAU) == 1)
message(sprintf("Threshold (tau) = %.4f", TAU))

# ── Constraint table → INTOLERANT set ────────────────────────────────────────
constr_raw <- read_tsv(CONSTRAINT_PATH, show_col_types = FALSE)

if ("mane_select" %in% names(constr_raw)) {
  mane_only <- constr_raw %>% filter(mane_select == TRUE)
  if (nrow(mane_only) > 0) constr_raw <- mane_only
}

constr <- constr_raw %>%
  rename(.gene = gene, .loeuf = `lof.oe_ci.upper`, .pli = `lof.pLI`) %>%
  filter(!is.na(.gene), !is.na(.loeuf), !is.na(.pli)) %>%
  group_by(.gene) %>%
  summarise(loeuf = min(.loeuf), pli = max(.pli), .groups = "drop")

intolerant_genes <- constr %>%
  filter(loeuf < LOEUF_MAX, pli >= PLI_MIN) %>% pull(.gene)

# ── AD gene list ─────────────────────────────────────────────────────────────
ad_df <- read_csv(AD_GENES_PATH, show_col_types = FALSE)
ad_col <- intersect(c("hgnc_symbol", "gene", "Gene", "GENE",
                      "symbol", "Symbol", "x", "X"), names(ad_df))[1]
ad_genes <- unique(trimws(as.character(ad_df[[ad_col]])))
ad_genes <- ad_genes[!is.na(ad_genes) & ad_genes != ""]

# ── INTOLERANT ∩ AD ──────────────────────────────────────────────────────────
keep_genes <- intersect(intolerant_genes, ad_genes)
message(sprintf("Intolerant: %d   AD: %d   Intersection: %d",
                length(intolerant_genes), length(ad_genes), length(keep_genes)))

both <- both_all %>% filter(hgnc_symbol %in% keep_genes)

gnom <- both %>% filter(group == "gnomAD")
plp  <- both %>% filter(group == "ClinVar P/LP")
vus  <- both %>% filter(group == "ClinVar VUS")
greg <- both %>% filter(group == "GREGoR")

for (ds in DATASETS) {
  d <- both %>% filter(group == ds)
  if (nrow(d) > 0) {
    message(sprintf("  %-13s n=%5d  escape=%.1f%%  genes=%d",
                    ds, nrow(d), 100 * mean(d$predicted_label == "escape"),
                    n_distinct(d$hgnc_symbol)))
  } else {
    message(sprintf("  %-13s n=0", ds))
  }
}

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
# PANEL A — overlapping density histograms
# ══════════════════════════════════════════════════════════════════════════════
pA <- ggplot(both, aes(x = escape_probability, fill = group, color = group)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, alpha = 0.45,
                 position = "identity", linewidth = 0.3) +
  geom_vline(xintercept = TAU, linetype = "dashed",
             color = "grey25", linewidth = 1.0) +
  annotate("text", x = TAU + 0.02, y = Inf, vjust = 4,
           label = sprintf("tau = %.2f", TAU),
           fontface = "bold.italic", size = 4.3, hjust = 0, color = "grey20") +
  scale_fill_manual(values = COL) +
  scale_color_manual(values = COL) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Escape probability", y = "Density", tag = "A") +
  fig_theme + theme(legend.position = c(0.88, 0.78))

# ══════════════════════════════════════════════════════════════════════════════
# PANEL B — class composition
# ══════════════════════════════════════════════════════════════════════════════
panel_b_df <- both %>%
  count(group, predicted_label) %>%
  group_by(group) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  mutate(class = if_else(predicted_label == "NMD", "NMD-sens.", "Escape"),
         class = factor(class, levels = c("NMD-sens.", "Escape")))

B_YMAX <- min(100, ceiling(max(panel_b_df$pct) / 10) * 10 + 10)

pB <- ggplot(panel_b_df, aes(x = class, y = pct, fill = group)) +
  geom_col(position = position_dodge(0.8), width = 0.75,
           color = "grey15", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", pct)),
            position = position_dodge(0.8), vjust = -0.4,
            fontface = "bold", size = 3.4) +
  scale_fill_manual(values = COL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), limits = c(0, B_YMAX)) +
  labs(x = NULL, y = "Variants (%)", tag = "B") +
  fig_theme + theme(legend.position = "top", legend.justification = "right")

# ══════════════════════════════════════════════════════════════════════════════
# PANEL C — escape fraction by chromosome
# ══════════════════════════════════════════════════════════════════════════════
chrom_levels <- paste0("chr", 1:22)
panel_c_df <- both %>%
  filter(CHROM %in% chrom_levels) %>%
  mutate(CHROM = factor(CHROM, levels = chrom_levels)) %>%
  group_by(group, CHROM) %>%
  summarise(pct_escape = 100 * mean(predicted_label == "escape"),
            .groups = "drop")

mean_lines <- both %>% group_by(group) %>%
  summarise(mean_pct = 100 * mean(predicted_label == "escape"), .groups = "drop")

pC <- ggplot(panel_c_df, aes(x = CHROM, y = pct_escape, fill = group)) +
  geom_col(position = position_dodge(0.85), width = 0.78,
           color = "grey15", linewidth = 0.3) +
  geom_hline(data = mean_lines, aes(yintercept = mean_pct, color = group),
             linetype = "dotted", linewidth = 0.9, show.legend = FALSE) +
  scale_fill_manual(values = COL) +
  scale_color_manual(values = COL) +
  scale_x_discrete(labels = 1:22, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Chromosome", y = "Escape calls (%)", tag = "C") +
  fig_theme + theme(legend.position = "top", legend.justification = "left",
                    legend.direction = "horizontal")   # above bars (no title now)

# ══════════════════════════════════════════════════════════════════════════════
# PANELS D-G — variant density per gene = variants / CDS length (per kb)
#   Stacked NMD-sensitive / Escape. Top 15 genes ranked by density among genes
#   with at least `min_n` variants (as in Figure 7).
# ══════════════════════════════════════════════════════════════════════════════

# ── CDS length per gene (same file as Figure 7) ──────────────────────────────
if (!file.exists(CDS_FILE)) {
  library(biomaRt)
  genes <- unique(na.omit(both_all$hgnc_symbol))
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  tx <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id",
                             "transcript_biotype", "transcript_is_canonical"),
              filters = "hgnc_symbol", values = genes, mart = mart) %>%
    filter(transcript_biotype == "protein_coding")
  cds <- getBM(attributes = c("ensembl_transcript_id", "cds_length"),
               filters = "ensembl_transcript_id",
               values = unique(tx$ensembl_transcript_id), mart = mart) %>%
    distinct() %>% filter(!is.na(cds_length))
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
                                min_n = MIN_N_DEFAULT) {
  df <- df %>% filter(!is.na(hgnc_symbol), hgnc_symbol != "")
  if (nrow(df) == 0)
    return(ggplot() +
             labs(tag = tag) + fig_theme)

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

  # the AD-intolerant subset is small: lower min_n until 15 genes qualify (>= 1)
  while (min_n > 1 && sum(gene_tot$total >= min_n) < 15) min_n <- min_n - 1
  message(sprintf("%s: min_n = %d (%d genes qualify)", dataset_name, min_n,
                  sum(gene_tot$total >= min_n)))

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
                                 "Escape"        = unname(col_light)),
                      drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
    labs(x = "Variants per kb of CDS", y = NULL, tag = tag) +
    fig_theme +
    theme(axis.text.y     = element_text(face = "bold.italic", size = FONT),
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
                          min_n = 3)   # GREGoR is small (as in Figure 7)

# ══════════════════════════════════════════════════════════════════════════════
# Combine
# ══════════════════════════════════════════════════════════════════════════════
combined <- (pA | pB) /
             pC        /
            (pD | pE) /
            (pF | pG) +
  plot_layout(heights = c(1.0, 0.9, 1.4, 1.4))

# cairo_pdf so the ≥ / — characters render in the PDF
ggsave("NMD_escape_AD_intolerant_panels_A_to_G_corrected.pdf", combined,
       width = 17, height = 22, device = cairo_pdf)
ggsave("NMD_escape_AD_intolerant_panels_A_to_G_corrected.png", combined,
       width = 17, height = 22, dpi = 300)

cat("\nSaved: NMD_escape_AD_intolerant_panels_A_to_G_corrected.pdf / .png\n")
print(combined)
