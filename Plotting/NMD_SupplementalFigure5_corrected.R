# ══════════════════════════════════════════════════════════════════════════════
# NMD — Supplemental Figure: 3'UTR Features vs NMD Efficiency
#
#   A = 3'UTR length (all variants) → Beeswarm + box, Short vs Long
#   B = 3'UTR length (last-exon)   → Beeswarm + box, Short vs Long, LAST-EXON ONLY
#   C = 3'UTR AU content           → Beeswarm + box by tertile
#   D = 3'UTR UC content           → Beeswarm + box by tertile
#   E = 3'UTR introns (last-exon)  → Violin + box (With vs Without), LAST-EXON ONLY
#   F = PhastCons conservation     → Bar ± CI by 3 categories
#   G = 3'UTR RBP motifs           → Lollipop: Bonferroni-sig only
#
# Panels A, C, D, F, G use the full variant set. Panels B and E use last-exon only.
# (GC content panel removed in this version.)
# Last-exon defined as last.EJC == "last.exon".
#
# Required file: TOPMed_stopgain_September25_corrected_readyformodel.csv
#
# v2 (corrected dataset, October 2026):
#   - input switched to the corrected September 25 file (as Figures 3-6, S2-S4)
#   - filter ALLELE.RAT >= 0.35 (was > 0.35), matching all other figures
#   - column check at load time; if the PhastCons category or log2 3'UTR length
#     columns were dropped in feature cleaning, they are rebuilt from the
#     numeric columns with the same cut-offs
#   - unused libraries (ggdist, gghalves) no longer loaded
#   - prints n per group for every panel (for the legend)
#   - output: SupplementalFigure5_revised.pdf / .png
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
# ggdist / gghalves were loaded but never used; removed in v2
library(ggbeeswarm)
library(patchwork)
library(scales)

select    <- dplyr::select; filter    <- dplyr::filter
mutate    <- dplyr::mutate; summarise <- dplyr::summarise
group_by  <- dplyr::group_by; distinct  <- dplyr::distinct
left_join <- dplyr::left_join; arrange   <- dplyr::arrange

# (No bootstrap — deterministic per-transcript median used throughout.)
TITLE_COL <- "#2E6DA4"; FONT <- 72   # bumped from 56 — much larger for poster-grade

base_theme <- theme_classic(base_size = FONT) + theme(
  legend.position    = "none",
  axis.text.x        = element_text(size = FONT + 2, face = "bold",
                                     color = "grey10", lineheight = 1.1),
  axis.text.y        = element_text(size = FONT + 2, face = "bold",
                                     color = "grey15"),
  axis.title.x       = element_text(size = FONT + 5, face = "bold",
                                     color = "grey10", margin = margin(t = 14)),
  axis.title.y       = element_text(size = FONT + 5, face = "bold",
                                     color = "grey10", margin = margin(r = 14)),
  axis.line          = element_line(color = "grey15", linewidth = 2.2),   # bumped
  axis.ticks         = element_line(color = "grey25", linewidth = 1.8),   # bumped
  axis.ticks.length  = unit(0.35, "cm"),
  panel.grid.major   = element_line(color = "grey90", linewidth = 0.85),
  panel.grid.minor   = element_blank(),
  plot.background    = element_rect(fill = "white", color = NA),
  panel.background   = element_rect(fill = "white", color = NA),
  strip.text         = element_text(size = FONT + 2, face = "bold",
                                     color = "grey10", margin = margin(8, 6, 8, 6)),
  strip.background   = element_rect(fill = "grey94", color = "grey40",
                                     linewidth = 1.7),                    # bumped
  plot.tag           = element_text(size = FONT + 22, face = "bold",
                                     color = "grey10"),
  plot.title         = element_text(size = FONT + 5,  face = "bold",
                                     color = TITLE_COL, hjust = 0.5,
                                     margin = margin(b = 12)),
  plot.caption       = element_text(size = FONT - 4,  face = "bold.italic",
                                     color = "grey50",
                                     hjust = 0.5, margin = margin(t = 10)),
  plot.margin        = margin(18, 22, 18, 22),
  panel.spacing      = unit(1.4, "lines"))

ref05      <- geom_hline(yintercept = 0.5, linetype = "dashed",
                          color = "grey30", linewidth = 2.0)               # bumped
y_nmd_full <- scale_y_continuous(breaks = c(0, .25, .5, .75, 1),
                                   limits = c(-.04, 1.58), expand = c(0, 0))
y_nmd_bar  <- scale_y_continuous(breaks = seq(0, 1, .25),
                                   limits = c(0, 1.05), expand = c(0, 0))

# Balanced quintile helper
make_q5 <- function(x) {
  cuts <- quantile(x, probs = seq(0, 1, .2), na.rm = TRUE)
  factor(paste0("Q", cut(x, breaks = cuts, include.lowest = TRUE,
                           labels = FALSE)),
         levels = paste0("Q", 1:5))
}

# Tertile helper
make_tertiles <- function(x) {
  cuts <- quantile(x, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  factor(c("Low","Medium","High")[cut(x, breaks = cuts,
                                       include.lowest = TRUE, labels = FALSE)],
         levels = c("Low","Medium","High"))
}

seq_comp <- function(lvls) lapply(seq_len(length(lvls) - 1),
                                   function(i) c(lvls[i], lvls[i + 1]))

# ══════════════════════════════════════════════════════════════════════════════
# LOAD & PREP
# ══════════════════════════════════════════════════════════════════════════════
INPUT_CSV <- "TOPMed_stopgain_September25_corrected_readyformodel.csv"
df_raw <- read.csv(INPUT_CSV, stringsAsFactors = FALSE) %>%
  filter(!is.na(ALLELE.RAT), !is.na(TxName),
         ALLELE.RAT >= 0.35)          # v2: >= (was >), same as all other figures

# ── v2: column check / rebuild from numeric columns if needed ───────────────
if (!"phastcons_new3utr_first200_median_cat" %in% names(df_raw)) {
  if (!"phastcons_new3utr_first200_median" %in% names(df_raw))
    stop("No 3'UTR PhastCons column found (expected phastcons_new3utr_first200_median[_cat]).")
  cat("Rebuilding phastcons_new3utr_first200_median_cat from the numeric column\n")
  df_raw$phastcons_new3utr_first200_median_cat <- with(df_raw, case_when(
    phastcons_new3utr_first200_median <  0.1 ~ "Low (<0.1)",
    phastcons_new3utr_first200_median <  0.5 ~ "Medium (0.1=< & <0.5)",
    phastcons_new3utr_first200_median >= 0.5 ~ "High (>=0.5)"))
}
if (!"log2_3utr" %in% names(df_raw) && "threeUTR_length" %in% names(df_raw)) {
  cat("Rebuilding log2_3utr = log2(threeUTR_length)\n")
  df_raw$log2_3utr <- log2(df_raw$threeUTR_length)
}
if (!"Freq.cat" %in% names(df_raw)) df_raw$Freq.cat <- NA_character_   # not plotted
req3 <- c("threeUTR_length", "log2_3utr", "threeUTR_AU_content",
          "threeUTR_UC_content", "threeUTR.introns", "last.EJC")
miss3 <- setdiff(req3, names(df_raw))
if (length(miss3) > 0)
  stop("Missing column(s): ", paste(miss3, collapse = ", "), "\nPossible matches: ",
       paste(grep("3utr|threeutr|utr3|last", names(df_raw), ignore.case = TRUE,
                  value = TRUE), collapse = ", "))
cat("threeUTR.introns values:\n"); print(table(df_raw$threeUTR.introns, useNA = "ifany"))

df_raw <- df_raw %>%
  mutate(
    # PhastCons conservation category (new UTR, first 200nt)
    phast_utr3 = case_when(
      phastcons_new3utr_first200_median_cat == "Low (<0.1)"             ~ "Low\n(<0.1)",
      phastcons_new3utr_first200_median_cat == "Medium (0.1=< & <0.5)" ~ "Medium\n(0.1\u20130.5)",
      phastcons_new3utr_first200_median_cat == "High (>=0.5)"          ~ "High\n(\u22650.5)",
      TRUE ~ NA_character_
    ),
    phast_utr3 = factor(phast_utr3,
                         levels = c("Low\n(<0.1)", "Medium\n(0.1\u20130.5)",
                                    "High\n(\u22650.5)"))
  )

utr3_cols <- grep("^utr3_all\\.", names(df_raw), value = TRUE)
cat(sprintf("Total utr3_all RBP cols: %d\n", length(utr3_cols)))
cat(sprintf("[Filter applied: ALLELE.RAT >= 0.35]\n"))
cat(sprintf("Variants after filter: %d  |  Unique tx: %d\n",
            nrow(df_raw), n_distinct(df_raw$TxName)))

# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
# TRANSCRIPT-LEVEL SUMMARY (deterministic — replaces prior bootstrap)
# For each transcript: median ALLELE.RAT across its surviving variants.
#   - 1 variant   → use that single value
#   - >1 variants → use the median across them
# ══════════════════════════════════════════════════════════════════════════════
cat("Computing transcript-level medians...\n")

tx_summary <- df_raw %>%
  group_by(TxName) %>%
  summarise(
    n_variants = n(),
    median_NMD = median(ALLELE.RAT, na.rm = TRUE),
    .groups    = "drop"
  )

cat(sprintf("  Total transcripts: %d  |  with >1 variant: %d (%.1f%%)\n",
            nrow(tx_summary),
            sum(tx_summary$n_variants > 1),
            100 * mean(tx_summary$n_variants > 1)))

tx_feat <- df_raw %>% distinct(TxName, .keep_all = TRUE) %>%
  mutate(
    # Derive utr_intron here so it survives the join correctly
    utr_intron = if_else(
      threeUTR.introns == "There is a 3UTR intron",
      "With 3\u2019UTR intron",
      "Without 3\u2019UTR intron"
    ),
    utr_intron = factor(utr_intron,
                         levels = c("Without 3\u2019UTR intron",
                                    "With 3\u2019UTR intron")),
    freq_grp = case_when(
      Freq.cat == "Ultra-rare variants"  ~ "Ultra-rare",
      Freq.cat == "Rare/Common variants" ~ "Rare/Common",
      TRUE ~ NA_character_
    ),
    freq_grp = factor(freq_grp, levels = c("Ultra-rare","Rare/Common"))
  ) %>%
  dplyr::select(TxName, threeUTR_length, log2_3utr,
                threeUTR_AU_content, threeUTR_UC_content,
                utr_intron, phast_utr3, freq_grp)

tx_df <- tx_summary %>% left_join(tx_feat, by = "TxName") %>% drop_na(median_NMD)
cat(sprintf("Transcript-level table: %d transcripts\n\n", nrow(tx_df)))

# ══════════════════════════════════════════════════════════════════════════════
# LAST-EXON TRANSCRIPT-LEVEL SUMMARY (deterministic — replaces bootstrap)
# Used by Panels B (last-exon 3'UTR length) and F (last-exon introns).
# Last-exon status is a variant-level property; we filter first, then collapse
# to transcript-level via median.
# ══════════════════════════════════════════════════════════════════════════════
if (!"last.EJC" %in% names(df_raw)) {
  stop("Column 'last.EJC' not found. Last-exon panels cannot be built.")
}

df_lastexon <- df_raw %>% filter(last.EJC == "last.exon")
cat(sprintf("Last-exon variants = %d (%.1f%% of dataset)\n",
            nrow(df_lastexon),
            100 * nrow(df_lastexon) / nrow(df_raw)))

tx_summary_le <- df_lastexon %>%
  group_by(TxName) %>%
  summarise(
    n_variants = n(),
    median_NMD = median(ALLELE.RAT, na.rm = TRUE),
    .groups    = "drop"
  )

# Join all transcript-level features we'll need (log2_3utr for panel B,
# utr_intron for panel F). These are transcript-level properties so the join
# is safe regardless of which variant was sampled per transcript.
tx_df_le <- tx_summary_le %>%
  left_join(tx_feat %>% dplyr::select(TxName, log2_3utr, threeUTR_length,
                                      utr_intron),
            by = "TxName") %>%
  drop_na(median_NMD)

cat(sprintf("Last-exon transcript-level table: %d transcripts\n\n",
            nrow(tx_df_le)))

# ══════════════════════════════════════════════════════════════════════════════
# HELPER — BEESWARM + BOX content tertile panel (used for panels C / D / E)
# ══════════════════════════════════════════════════════════════════════════════
make_content_panel <- function(tx_df, col, x_label, title, tag,
                                pal_low, pal_med, pal_high) {
  df_p <- tx_df %>% drop_na(!!sym(col)) %>%
    mutate(grp = make_tertiles(.data[[col]]))

  pal3 <- c("Low" = pal_low, "Medium" = pal_med, "High" = pal_high)

  stat_p <- df_p %>%
    wilcox_test(median_NMD ~ grp,
                comparisons = list(c("Low","Medium"), c("Medium","High"),
                                   c("Low","High"))) %>%
    add_significance("p") %>% add_xy_position(x = "grp") %>%
    mutate(y.position = c(1.10, 1.20, 1.30), grp = group1)   # v2: inherited aes

  n_p <- df_p %>% count(grp) %>% mutate(label = format(n, big.mark = ","))

  ggplot(df_p, aes(x = grp, y = median_NMD, fill = grp, color = grp)) +
    geom_quasirandom(aes(color = grp), size = 1.4, alpha = 0.28,        # was 0.90 / 0.22
                     shape = 16, width = 0.28, bandwidth = 0.10) +
    geom_boxplot(aes(fill = grp), width = 0.24, outlier.shape = NA,
                 alpha = 0.55, linewidth = 2.2, color = "grey10",
                 fatten = 3.5) +                                          # was 2.5
    stat_summary(fun = median, geom = "point", size = 9, shape = 21,     # was 6
                 fill = "white", color = "black", stroke = 2.2) +        # stroke was 1.4
    ref05 +
    geom_text(data = n_p, aes(x = grp, y = 1.48, label = label),
              inherit.aes = FALSE, size = 14, fontface = "bold",         # bumped
              color = "grey15") +
    stat_pvalue_manual(stat_p, label = "p.signif",
                       tip.length = 0.012, bracket.size = 1.2,           # was 0.70
                       size = 18, color = "grey15") +                    # bumped
    scale_fill_manual(values = pal3) +
    scale_color_manual(values = pal3) +
    y_nmd_full +
    labs(x = x_label, y = "NMD efficiency",
         title = title, tag = tag) +
    base_theme
}

# ══════════════════════════════════════════════════════════════════════════════
# PANELS C, D — BEESWARM + BOX: AU / UC content by tertile
#               (GC content panel removed)
# ══════════════════════════════════════════════════════════════════════════════
pB <- make_content_panel(tx_df, "threeUTR_AU_content",
                          "3\u2019UTR AU content tertile",
                          "3\u2019UTR AU Content", "C",
                          "#D6EAF8","#2E86C1","#1A5276")

pD <- make_content_panel(tx_df, "threeUTR_UC_content",
                          "3\u2019UTR UC content tertile",
                          "3\u2019UTR UC Content", "D",
                          "#FDEBD0","#E67E22","#784212")

# ══════════════════════════════════════════════════════════════════════════════
# PANEL F — VIOLIN + BOX: 3'UTR introns (With vs Without), LAST-EXON only
# ══════════════════════════════════════════════════════════════════════════════
dfE <- tx_df_le %>% drop_na(utr_intron)

stat_E <- dfE %>%
  wilcox_test(median_NMD ~ utr_intron) %>%
  add_significance("p") %>% add_xy_position(x = "utr_intron") %>%
  mutate(y.position = 1.18, utr_intron = group1)
nE <- dfE %>% count(utr_intron) %>% mutate(label = format(n, big.mark = ","))

# v2: names now use the same curly apostrophe (\u2019) as the group labels;
# with a straight ' the colours never matched and panel E was drawn grey
pal2_E <- c("Without 3\u2019UTR intron" = "#AED6F1",
             "With 3\u2019UTR intron"    = "#C0392B")

pE <- ggplot(dfE, aes(x = utr_intron, y = median_NMD,
                       fill = utr_intron, color = utr_intron)) +
  geom_violin(alpha = 0.42, width = 0.80, trim = TRUE, linewidth = 1.8) +    # was 0.35
  geom_boxplot(aes(fill = utr_intron), width = 0.24, outlier.shape = 21,
               outlier.size = 2.2, outlier.alpha = 0.30,
               alpha = 0.55, linewidth = 2.2, color = "grey10",
               fatten = 3.5) +                                                # was 2.5
  stat_summary(fun = median, geom = "point", size = 9, shape = 21,           # was 6
               fill = "white", color = "black", stroke = 2.2) +              # stroke was 1.4
  ref05 +
  geom_text(data = nE, aes(x = utr_intron, y = 1.48, label = label),
            inherit.aes = FALSE, size = 14, fontface = "bold",               # bumped
            color = "grey15") +
  stat_pvalue_manual(stat_E, label = "p.signif",
                     tip.length = 0.014, bracket.size = 1.2,                 # was 0.70
                     size = 19, color = "grey15") +                          # bumped
  scale_fill_manual(values = pal2_E) +
  scale_color_manual(values = pal2_E) +
  y_nmd_full +
  labs(x = "3\u2019UTR intron status",
       y = "NMD efficiency",
       title = "3\u2019UTR Introns (Last-Exon Variants)",
       tag = "E") +
  base_theme

# ══════════════════════════════════════════════════════════════════════════════
# PANEL G — BAR ± CI: PhastCons conservation (3'UTR first 200 nt)
# ══════════════════════════════════════════════════════════════════════════════
dfF   <- tx_df %>% drop_na(phast_utr3)
bar_F <- dfF %>% group_by(phast_utr3) %>%
  summarise(mean_NMD = mean(median_NMD, na.rm = TRUE),
            se       = sd(median_NMD, na.rm = TRUE) / sqrt(n()),
            ci_lo    = mean_NMD - 1.96 * se,
            ci_hi    = mean_NMD + 1.96 * se,
            n        = n(), .groups = "drop")
nF    <- dfF %>% count(phast_utr3) %>% mutate(label = format(n, big.mark = ","))

stat_F <- dfF %>%
  wilcox_test(median_NMD ~ phast_utr3,
              comparisons = list(
                c("Low\n(<0.1)", "Medium\n(0.1\u20130.5)"),
                c("Medium\n(0.1\u20130.5)", "High\n(\u22650.5)"),
                c("Low\n(<0.1)", "High\n(\u22650.5)"))) %>%
  add_significance("p") %>% add_xy_position(x = "phast_utr3") %>%
  mutate(y.position = c(0.82, 0.90, 0.98), phast_utr3 = group1)

pal3_F <- c("Low\n(<0.1)"         = "#BDC3C7",
             "Medium\n(0.1\u20130.5)" = "#85C1E9",
             "High\n(\u22650.5)"  = "#1B4F72")

pF <- ggplot(bar_F, aes(x = phast_utr3, y = mean_NMD, fill = phast_utr3)) +
  geom_col(alpha = 0.88, color = "grey15", linewidth = 1.8, width = 0.72) +   # was 0.55
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.28, linewidth = 2.6, color = "grey10") +            # was 1.2
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "grey30", linewidth = 1.4) +                             # was 0.75
  geom_text(aes(label = sprintf("%.3f", mean_NMD), y = ci_hi + 0.020),
            size = 14, fontface = "bold", color = "grey10") +                 # bumped
  geom_text(data = nF, aes(x = phast_utr3, y = 0.04, label = label),
            inherit.aes = FALSE, size = 13, fontface = "bold",                # was 7.5
            color = "white") +
  stat_pvalue_manual(stat_F, label = "p.signif",
                     tip.length = 0.014, bracket.size = 1.2,                  # was 0.70
                     size = 18, color = "grey15") +                           # bumped
  scale_fill_manual(values = pal3_F) +
  y_nmd_bar +
  labs(x = "PhastCons conservation (3\u2019UTR first 200 nt)",
       y = "Mean NMD efficiency",
       title = "3\u2019UTR PhastCons Conservation",
       tag = "F") +
  base_theme

# ══════════════════════════════════════════════════════════════════════════════
# PANEL H — LOLLIPOP: Bonferroni-significant 3'UTR RBP motifs
#           TRANSCRIPT-LEVEL: each transcript contributes a single row
#           (median ALLELE.RAT) with its 3'UTR RBP-motif annotation, since
#           motif presence is determined by 3'UTR sequence and is constant
#           across all PTC variants on the same transcript.
# ══════════════════════════════════════════════════════════════════════════════
cat("Computing transcript-level 3'UTR RBP motif statistics...\n")

# Build transcript-level RBP table: median NMD per transcript joined with the
# (constant per transcript) RBP-motif flags taken from one row each.
tx_rbp <- tx_summary %>%
  left_join(
    df_raw %>% distinct(TxName, .keep_all = TRUE) %>%
      dplyr::select(TxName, all_of(utr3_cols)),
    by = "TxName"
  ) %>%
  drop_na(median_NMD)

cat(sprintf("  Transcript-level RBP table: %d transcripts\n", nrow(tx_rbp)))

utr3_stats <- map_dfr(utr3_cols, function(col) {
  has <- tx_rbp %>% filter(.data[[col]] == 1) %>% pull(median_NMD)
  no  <- tx_rbp %>% filter(.data[[col]] == 0) %>% pull(median_NMD)
  if (length(has) < 10) return(NULL)
  wt <- wilcox.test(has, no, alternative = "two.sided")
  tibble(rbp   = sub("utr3_all\\.", "", col),
         p_raw = wt$p.value,
         n_has = length(has),
         diff  = mean(has, na.rm = TRUE) - mean(no, na.rm = TRUE))
}) %>%
  mutate(
    p_bonf = pmin(p_raw * n(), 1),
    sig    = case_when(p_bonf < .0001 ~ "****", p_bonf < .001 ~ "***",
                        p_bonf < .01   ~ "**",  p_bonf < .05  ~ "*",
                        TRUE ~ "ns")
  )

sig_utr3 <- utr3_stats %>%
  filter(p_bonf < 0.05) %>%
  arrange(diff) %>%
  mutate(
    direction = if_else(diff < 0, "Inhibits NMD", "Promotes NMD"),
    rbp       = factor(rbp, levels = rbp)
  )

n_neg <- sum(sig_utr3$diff < 0)
n_pos <- sum(sig_utr3$diff > 0)
cat(sprintf("Bonferroni-significant 3'UTR RBPs (transcript-level): %d (inhibiting=%d, promoting=%d)\n",
            nrow(sig_utr3), n_neg, n_pos))
print(sig_utr3 %>% dplyr::select(rbp, p_bonf, n_has, diff, sig))

x_range <- range(sig_utr3$diff)
# Always include zero in the visible range so segments anchor at the y-axis
# (otherwise, when all RBPs are on one side, segments extend off-panel and
#  only the endpoint dots appear visible).
x_lo_raw <- min(0, x_range[1])
x_hi_raw <- max(0, x_range[2])
x_pad   <- (x_hi_raw - x_lo_raw) * 0.30
x_lo    <- x_lo_raw - x_pad
x_hi    <- x_hi_raw + x_pad
n_rows  <- nrow(sig_utr3)

pG <- ggplot(sig_utr3, aes(y = rbp)) +

  geom_vline(xintercept = 0, color = "grey20", linewidth = 1.6) +              # was 0.9
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#EBF5FB", alpha = 0.45) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#EAFAF1", alpha = 0.45) +

  # Lollipop segment: a thin dark underline first for definition,
  # then the colored segment on top, fully opaque.
  geom_segment(aes(x = 0, xend = diff, y = rbp, yend = rbp),
               linewidth = 9.5, color = "grey15", alpha = 1.0) +
  geom_segment(aes(x = 0, xend = diff, y = rbp, yend = rbp, color = direction),
               linewidth = 7.0, alpha = 1.0) +
  geom_point(aes(x = diff, color = direction), size = 14, alpha = 1.0,
             stroke = 1.5) +

  geom_text(aes(x = diff,
                label = paste0("n=", format(n_has, big.mark = ",")),
                hjust = if_else(diff < 0, 1.15, -0.15)),
            size = 13, fontface = "italic", color = "grey30") +                # bumped
  geom_text(aes(x = diff + sign(diff) * (x_hi_raw - x_lo_raw) * 0.012,
                label = sig,
                hjust = if_else(diff < 0, -0.15, 1.15)),
            size = 17, fontface = "bold", color = "grey10") +                  # bumped

  {if (n_neg > 0) annotate("text", x = x_lo_raw / 2, y = n_rows + 1.2,
                             label = "Inhibits NMD", fontface = "bold.italic",
                             size = 17, color = "#2980B9", hjust = 0.5)} +     # bumped
  {if (n_pos > 0) annotate("text", x = x_hi_raw / 2, y = n_rows + 1.2,
                             label = "Promotes NMD", fontface = "bold.italic",
                             size = 17, color = "#27AE60", hjust = 0.5)} +     # bumped

  scale_color_manual(values = c("Inhibits NMD" = "#2980B9",
                                 "Promotes NMD" = "#27AE60"),
                     guide = "none") +
  scale_x_continuous(breaks = pretty(c(x_lo, x_hi), n = 6),
                     limits = c(x_lo, x_hi),
                     labels = function(x) sprintf("%.2f", x),
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off") +

  labs(x       = "Mean NMD efficiency difference (motif present \u2212 absent)",
       y       = NULL,
       title   = "3\u2019UTR RBP Binding Motifs vs NMD Efficiency",
       tag     = "G") +
  base_theme +
  theme(
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.major.x = element_line(color = "grey91", linewidth = 0.45),
    axis.text.y        = element_text(size = FONT + 1, face = "bold",
                                       color = "grey10"),
    plot.margin        = margin(24, 36, 16, 20)
  )

# ══════════════════════════════════════════════════════════════════════════════
# PANEL A — BEESWARM + BOX: 3'UTR length (Short/Long), all variants
# ══════════════════════════════════════════════════════════════════════════════
# Fixed log2 threshold = 10  (i.e. 2^10 = 1,024 nt)
UTR_LOG2_THRESH <- 10

cat(sprintf("3\'UTR length cut: log2 = %.0f  (= %.0f nt)\n",
            UTR_LOG2_THRESH, 2^UTR_LOG2_THRESH))

dfA <- tx_df %>% drop_na(log2_3utr) %>%
  mutate(grp = factor(
    if_else(log2_3utr < UTR_LOG2_THRESH,
            "Short
(log₂ 3\'UTR < 10)",
            "Long
(log₂ 3\'UTR ≥ 10)"),
    levels = c("Short
(log₂ 3\'UTR < 10)",
               "Long
(log₂ 3\'UTR ≥ 10)")))

pal2_A <- c("Short
(log₂ 3\'UTR < 10)" = "#C39BD3",
            "Long
(log₂ 3\'UTR ≥ 10)"  = "#4A235A")

nA <- dfA %>% count(grp) %>% mutate(label = format(n, big.mark = ","))

stat_A <- dfA %>%
  wilcox_test(median_NMD ~ grp) %>%
  add_significance("p") %>% add_xy_position(x = "grp") %>%
  mutate(y.position = 1.18, grp = group1)

pA <- ggplot(dfA, aes(x = grp, y = median_NMD, fill = grp, color = grp)) +
  geom_quasirandom(aes(color = grp), size = 1.4, alpha = 0.28,               # was 0.90 / 0.22
                   shape = 16, width = 0.28, bandwidth = 0.10) +
  geom_boxplot(aes(fill = grp), width = 0.24, outlier.shape = NA,
               alpha = 0.55, linewidth = 2.2, color = "grey10",
               fatten = 3.5) +                                                  # was 2.5
  stat_summary(fun = median, geom = "point", size = 9, shape = 21,            # was 6
               fill = "white", color = "black", stroke = 2.2) +               # stroke was 1.4
  ref05 +
  geom_text(data = nA, aes(x = grp, y = 1.48, label = label),
            inherit.aes = FALSE, size = 14, fontface = "bold",                 # bumped
            color = "grey15") +
  stat_pvalue_manual(stat_A, label = "p.signif",
                     tip.length = 0.012, bracket.size = 1.2,                   # was 0.70
                     size = 18, color = "grey15") +                            # bumped
  scale_fill_manual(values = pal2_A) +
  scale_color_manual(values = pal2_A) +
  y_nmd_full +
  labs(x = "3\u2019UTR length (log\u2082 scale)",
       y = "NMD efficiency",
       title = "3\u2019UTR Length",
       tag = "A") +
  base_theme

# ══════════════════════════════════════════════════════════════════════════════
# PANEL B — BEESWARM + BOX: 3'UTR length (Short/Long), LAST-EXON variants only
# ══════════════════════════════════════════════════════════════════════════════
# Same log2 threshold as panel A (log2 3'UTR = 10 -> 1,024 nt), but restricted
# to variants where last.EJC == "last.exon" (PTCs that largely evade NMD).
dfB_le <- tx_df_le %>% drop_na(log2_3utr) %>%
  mutate(grp = factor(
    if_else(log2_3utr < UTR_LOG2_THRESH,
            "Short\n(log\u2082 3\u2019UTR < 10)",
            "Long\n(log\u2082 3\u2019UTR \u2265 10)"),
    levels = c("Short\n(log\u2082 3\u2019UTR < 10)",
               "Long\n(log\u2082 3\u2019UTR \u2265 10)")))

pal2_B_le <- c("Short\n(log\u2082 3\u2019UTR < 10)" = "#F5B7B1",
               "Long\n(log\u2082 3\u2019UTR \u2265 10)" = "#7B241C")

nB_le <- dfB_le %>% count(grp) %>% mutate(label = format(n, big.mark = ","))

stat_B_le <- dfB_le %>%
  wilcox_test(median_NMD ~ grp) %>%
  add_significance("p") %>% add_xy_position(x = "grp") %>%
  mutate(y.position = 1.18, grp = group1)

pB_le <- ggplot(dfB_le, aes(x = grp, y = median_NMD, fill = grp, color = grp)) +
  geom_quasirandom(aes(color = grp), size = 1.4, alpha = 0.28,
                   shape = 16, width = 0.28, bandwidth = 0.10) +
  geom_boxplot(aes(fill = grp), width = 0.24, outlier.shape = NA,
               alpha = 0.55, linewidth = 2.2, color = "grey10",
               fatten = 3.5) +
  stat_summary(fun = median, geom = "point", size = 9, shape = 21,
               fill = "white", color = "black", stroke = 2.2) +
  ref05 +
  geom_text(data = nB_le, aes(x = grp, y = 1.48, label = label),
            inherit.aes = FALSE, size = 14, fontface = "bold",
            color = "grey15") +
  stat_pvalue_manual(stat_B_le, label = "p.signif",
                     tip.length = 0.014, bracket.size = 1.6,
                     size = 22, color = "grey15") +
  scale_fill_manual(values = pal2_B_le) +
  scale_color_manual(values = pal2_B_le) +
  y_nmd_full +
  labs(x = "3\u2019UTR length (log\u2082 scale)",
       y = "NMD efficiency",
       title = "3\u2019UTR Length (Last-Exon Variants)",
       tag = "B") +
  base_theme

# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
# COMBINE
#   Row 1: A (3'UTR length, all)        |  B (3'UTR length, last-exon)
#   Row 2: C (AU content)               |  D (UC content)
#   Row 3: E (introns)                  |  F (PhastCons conservation)
#   Row 4: G (RBP lollipop, full width)
# ══════════════════════════════════════════════════════════════════════════════
row1 <- (pA | pB_le) + plot_layout(widths = c(1.0, 1.0))

combined <- row1 / (pB | pD) / (pE | pF) / pG +
  plot_layout(heights = c(0.9, 1, 1, 0.7)) &
  theme(plot.margin = margin(12, 16, 12, 16))

combined <- combined +
  plot_annotation(
    title   = "3\u2019UTR Features vs NMD Efficiency",
    theme   = theme(
      plot.title   = element_text(size = FONT + 10, face = "bold",            # was +6
                                   color = TITLE_COL, hjust = 0.5,
                                   margin = margin(b = 10, t = 10)),
      plot.caption = element_text(size = FONT - 2, face = "bold.italic",       # was -4
                                   color = "grey50", hjust = 0.5,
                                   margin = margin(t = 12, b = 8))))

# ── v2: n per group for the legend ─────────────────────────────────────────
cat("\n── n per group (transcripts) ──\n")
show_n <- function(lbl, d, g) { cat(lbl, ": ", paste(paste0(gsub("\n", " ", names(table(d[[g]]))),
                                 " = ", table(d[[g]])), collapse = "; "), "\n", sep = "") }
show_n("A 3'UTR length (all)",       dfA,    "grp")
show_n("B 3'UTR length (last exon)", dfB_le, "grp")
show_n("E 3'UTR intron (last exon)", dfE,    "utr_intron")
show_n("F PhastCons",                dfF,    "phast_utr3")
cat("C/D: AU and UC tertiles, see plot labels\n")

ggsave("SupplementalFigure5_revised.pdf", combined,
       width = 60, height = 78, dpi = 300, limitsize = FALSE)
ggsave("SupplementalFigure5_revised.png", combined,
       width = 60, height = 78, dpi = 150, limitsize = FALSE)

cat("Saved: SupplementalFigure5_revised.pdf / .png\n")
