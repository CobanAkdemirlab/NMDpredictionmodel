# ══════════════════════════════════════════════════════════════════════════════
# NMD — Transcript-Level Features vs NMD Efficiency  (v5)
#
#   A = Exon count             → Violin + box by tertile
#   B = Median expression      → Violin + box by tertile (log2 TPM)
#   C = mRNA half-life PC1     → Violin + box by tertile
#   D = Isoform count          → Violin + box, 3 SHAP-aligned categories
#                                  (1-5 / 6-11 / >=12)
#   E = CDS AU content         → Scatter + GAM
#   F = CDS UC content         → Scatter + GAM
#
# Changes vs v4:
#   1) Removed the two panels with only non-significant comparisons
#      (Reviewer 2, major comment 2):
#        - CDS length quintiles (was A)
#        - Whole Blood expression tertiles (was D)
#   2) Panels A-D now share ONE plot style (violin + box + median point,
#      n labels, sequential Wilcoxon brackets), replacing the previous mix of
#      violin / beeswarm / bar / sina / plain boxplot. Scatter panels unchanged.
#   3) Re-lettered A-F; layout 3 x 2.
#   4) Input: September 25 data frame; filter ALLELE.RAT >= 0.35, matching
#      Figures 3, 4 and 6.
#   5) cairo_pdf output (renders the Greek rho reliably); p-value tables carry
#      the grouping column so stat_pvalue_manual() works on all ggplot2/ggpubr
#      versions.
#
# Required file: TOPMed_stopgain_September25_corrected_readyformodel.csv
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)
library(scales)

select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
distinct  <- dplyr::distinct
left_join <- dplyr::left_join

set.seed(42)

# ── PDF device that works on any machine ─────────────────────────────────────
#   macOS      -> built-in Quartz PDF (no XQuartz needed; renders rho, +/-)
#   elsewhere  -> cairo_pdf if available, otherwise pdf()
pdf_device <- function(filename, width, height, ...) {
  if (Sys.info()[["sysname"]] == "Darwin") {
    # macOS: the built-in Quartz PDF device needs no XQuartz. (cairo_pdf on a
    # Mac without XQuartz fails with "failed to load cairo DLL" even though
    # capabilities("cairo") can still report TRUE.)
    grDevices::quartz(type = "pdf", file = filename,
                      width = width, height = height, bg = "white")
  } else if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(filename = filename, width = width, height = height)
  } else {
    grDevices::pdf(file = filename, width = width, height = height)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# SHARED THEME (unchanged from v4)
# ══════════════════════════════════════════════════════════════════════════════
TITLE_COL <- "#2E6DA4"
FONT      <- 38

base_theme <- theme_classic(base_size = FONT) +
  theme(
    legend.position    = "none",
    axis.text.x        = element_text(size = FONT + 2, face = "bold",
                                      color = "grey10", lineheight = 1.1),
    axis.text.y        = element_text(size = FONT + 2, face = "bold",
                                      color = "grey15"),
    axis.title.x       = element_text(size = FONT + 5, face = "bold",
                                      color = "grey10", margin = margin(t = 14)),
    axis.title.y       = element_text(size = FONT + 5, face = "bold",
                                      color = "grey10", margin = margin(r = 14)),
    axis.line          = element_line(color = "grey15", linewidth = 1.8),
    axis.ticks         = element_line(color = "grey25", linewidth = 1.4),
    axis.ticks.length  = unit(0.38, "cm"),
    panel.grid.major   = element_line(color = "grey90", linewidth = 0.7),
    panel.grid.minor   = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.tag           = element_text(size = FONT + 28, face = "bold",
                                      color = "grey10"),
    plot.title         = element_text(size = FONT + 14, face = "bold",
                                      color = TITLE_COL, hjust = 0.5,
                                      margin = margin(b = 14)),
    plot.caption       = element_blank(),
    plot.subtitle      = element_blank(),
    plot.margin        = margin(18, 22, 18, 22),
    panel.spacing      = unit(1.4, "lines")
  )

scatter_theme <- base_theme +
  theme(panel.grid.major.x = element_line(color = "grey90", linewidth = 0.7))

ref05      <- geom_hline(yintercept = 0.5, linetype = "dashed",
                         color = "grey30", linewidth = 1.5)
y_nmd_full <- scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                                 limits = c(-0.04, 1.56), expand = c(0, 0))
y_nmd_scat <- scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                                 limits = c(0, 1.05), expand = c(0, 0))

# ══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════════════
make_tertiles <- function(x, prefix = "T", fmt = "%.0f") {
  cuts <- quantile(x, probs = seq(0, 1, 1/3), na.rm = TRUE)
  grp  <- cut(x, breaks = cuts, include.lowest = TRUE, labels = FALSE)
  labs <- sprintf(paste0(prefix, "%d\n(", fmt, "–", fmt, ")"),
                  1:3, cuts[1:3], cuts[2:4])
  factor(labs[grp], levels = labs)
}

seq_comparisons <- function(lvls)
  lapply(seq_len(length(lvls) - 1), function(i) c(lvls[i], lvls[i + 1]))

# One shared style for every categorical panel:
#   violin (light fill, dark outline) + box + white median point + n labels
#   + Wilcoxon brackets
make_violin_panel <- function(df, pal_dark, pal_light, comparisons,
                              x_lab, title, tag) {
  names(pal_dark)  <- levels(df$grp)
  names(pal_light) <- levels(df$grp)
  n_lab <- df %>% count(grp) %>% mutate(label = format(n, big.mark = ","))
  y_pos <- c(1.10, 1.20, 1.30)[seq_along(comparisons)]
  stat  <- df %>%
    wilcox_test(median_NMD ~ grp, comparisons = comparisons) %>%
    add_significance("p") %>% add_xy_position(x = "grp") %>%
    mutate(y.position = y_pos, grp = group1)

  ggplot(df, aes(x = grp, y = median_NMD, fill = grp, color = grp)) +
    geom_violin(alpha = 0.45, width = 0.80, trim = TRUE, linewidth = 1.2) +
    geom_boxplot(width = 0.24, outlier.shape = NA, alpha = 0.55,
                 linewidth = 1.3, color = "grey10", fatten = 3.5) +
    stat_summary(fun = median, geom = "point", size = 8, shape = 21,
                 fill = "white", color = "black", stroke = 2.0) +
    ref05 +
    geom_text(data = n_lab, aes(x = grp, y = 1.47, label = label),
              inherit.aes = FALSE, size = 11, fontface = "bold",
              color = "grey10") +
    stat_pvalue_manual(stat, label = "p.signif", tip.length = 0.012,
                       bracket.size = 1.1, size = 13, color = "grey15") +
    scale_fill_manual(values = pal_light) +
    scale_color_manual(values = pal_dark) +
    y_nmd_full +
    labs(x = x_lab, y = "NMD efficiency", title = title, tag = tag) +
    base_theme
}

# ══════════════════════════════════════════════════════════════════════════════
# LOAD & PREP — ALLELE.RAT >= 0.35; transcript level = median per transcript
# ══════════════════════════════════════════════════════════════════════════════
df_raw <- read.csv("TOPMed_stopgain_September25_corrected_readyformodel.csv",
                   stringsAsFactors = FALSE) %>%
  filter(!is.na(ALLELE.RAT), ALLELE.RAT >= 0.35, !is.na(TxName))
cat(sprintf("Variants after ALLELE.RAT >= 0.35 filter: %d\n", nrow(df_raw)))

tx_nmd <- df_raw %>%
  group_by(TxName) %>%
  summarise(median_NMD = median(ALLELE.RAT, na.rm = TRUE),
            n_var = dplyr::n(), .groups = "drop")

# Transcript-level features: first non-missing value per transcript
tx_features <- df_raw %>%
  group_by(TxName) %>%
  summarise(across(c(exon_count, MedianExpression_log2, half_life_PC1,
                     isoform_count, cdsseqs_AU_content, cdsseqs_UC_content),
                   ~ dplyr::first(na.omit(.x))),
            .groups = "drop")

tx_df <- tx_nmd %>% left_join(tx_features, by = "TxName") %>%
  drop_na(median_NMD)
cat(sprintf("Transcript-level table: %d transcripts\n\n", nrow(tx_df)))

# ══════════════════════════════════════════════════════════════════════════════
# PANEL A — Exon count tertiles
# ══════════════════════════════════════════════════════════════════════════════
dfA <- tx_df %>% drop_na(exon_count) %>%
  mutate(grp = make_tertiles(exon_count, "T", "%.0f"))
pA <- make_violin_panel(dfA,
                        pal_dark  = c("#6BB8D4", "#F0A500", "#D05050"),
                        pal_light = c("#A8D8EA", "#F8D07A", "#EEA0A0"),
                        comparisons = seq_comparisons(levels(dfA$grp)),
                        x_lab = "Exon count tertile",
                        title = "Exon Count", tag = "A")

# ══════════════════════════════════════════════════════════════════════════════
# PANEL B — Median expression tertiles
# ══════════════════════════════════════════════════════════════════════════════
dfB <- tx_df %>% drop_na(MedianExpression_log2) %>%
  mutate(grp = make_tertiles(MedianExpression_log2, "T", "%.1f"))
pB <- make_violin_panel(dfB,
                        pal_dark  = c("#4CAF7D", "#8E44AD", "#E67E22"),
                        pal_light = c("#A9DCC0", "#CDA6DD", "#F4C095"),
                        comparisons = seq_comparisons(levels(dfB$grp)),
                        x_lab = "Median expression tertile (log2 TPM)",
                        title = "Median Expression", tag = "B")

# ══════════════════════════════════════════════════════════════════════════════
# PANEL C — mRNA half-life PC1 tertiles
# ══════════════════════════════════════════════════════════════════════════════
dfC <- tx_df %>% drop_na(half_life_PC1) %>%
  mutate(grp = make_tertiles(half_life_PC1, "T", "%.1f"))
pC <- make_violin_panel(dfC,
                        pal_dark  = c("#E74C3C", "#F0A500", "#2980B9"),
                        pal_light = c("#F3A59D", "#F8D07A", "#A9CCE3"),
                        comparisons = seq_comparisons(levels(dfC$grp)),
                        x_lab = "mRNA half-life PC1 tertile",
                        title = "mRNA Half-Life PC1", tag = "C")

# ══════════════════════════════════════════════════════════════════════════════
# PANEL D — Isoform count, 3 SHAP-aligned categories
#   T1 = 1-5, T2 = 6-11, T3 = >=12 (SHAP zero-crossing ~13)
# ══════════════════════════════════════════════════════════════════════════════
lv_D <- c("T1\n(1–5)", "T2\n(6–11)", "T3\n(≥12)")
dfD <- tx_df %>% drop_na(isoform_count) %>%
  mutate(grp = factor(case_when(isoform_count <=  5 ~ lv_D[1],
                                isoform_count <= 11 ~ lv_D[2],
                                TRUE                ~ lv_D[3]),
                      levels = lv_D))
pD <- make_violin_panel(dfD,
                        pal_dark  = c("#A569BD", "#6C3483", "#4A235A"),
                        pal_light = c("#D7BDE2", "#B07CC6", "#8E5BA3"),
                        comparisons = list(lv_D[1:2], lv_D[2:3], lv_D[c(1, 3)]),
                        x_lab = "Isoform count category",
                        title = "Isoform Count", tag = "D")

# ══════════════════════════════════════════════════════════════════════════════
# PANELS E, F — Scatter + GAM (CDS AU, CDS UC) — style unchanged
# ══════════════════════════════════════════════════════════════════════════════
make_content_scatter <- function(data, x_col, x_label, title, tag, pt_col) {
  sub <- data %>%
    dplyr::select(x_val = !!sym(x_col), y_val = median_NMD) %>%
    drop_na()
  r_v   <- cor(sub$x_val, sub$y_val, use = "complete.obs", method = "spearman")
  x_rng <- range(sub$x_val)
  x_r   <- x_rng[1] + diff(x_rng) * 0.60
  x_n   <- x_rng[1] + diff(x_rng) * 0.02

  ggplot(sub, aes(x = x_val, y = y_val)) +
    geom_point(color = pt_col, size = 6.0, alpha = 0.42, shape = 16) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                color = "black", linewidth = 1.9,
                fill = "grey72", alpha = 0.25, se = TRUE) +
    ref05 +
    annotate("text", x = x_r, y = 0.06,
             label = sprintf("ρ = %.3f", r_v),
             size = 11, fontface = "bold.italic", color = "grey15") +
    annotate("text", x = x_n, y = 0.98,
             label = paste0("n = ", format(nrow(sub), big.mark = ",")),
             hjust = 0, size = 11, fontface = "bold", color = "grey15") +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       expand = c(0.03, 0)) +
    y_nmd_scat +
    labs(x = x_label, y = "NMD efficiency", title = title, tag = tag) +
    scatter_theme
}

pE <- make_content_scatter(tx_df, "cdsseqs_AU_content",
                           "CDS AU content", "CDS AU Content", "E", "#27AE60")
pF <- make_content_scatter(tx_df, "cdsseqs_UC_content",
                           "CDS UC content", "CDS UC Content", "F", "#7F8C8D")

# ══════════════════════════════════════════════════════════════════════════════
# COMBINE — 3 x 2
#   Row 1: A | B | C
#   Row 2: D | E | F
# ══════════════════════════════════════════════════════════════════════════════
combined <- (pA | pB | pC) / (pD | pE | pF) &
  theme(plot.margin = margin(12, 16, 12, 16))

ggsave("NMD_TranscriptFeatures_Figure_v5.pdf", combined, device = pdf_device,
       width = 46, height = 31, limitsize = FALSE)
ggsave("NMD_TranscriptFeatures_Figure_v5.png", combined,
       width = 46, height = 31, dpi = 200, limitsize = FALSE)
cat("Saved: NMD_TranscriptFeatures_Figure_v5.pdf / .png\n")
