# ══════════════════════════════════════════════════════════════════════════════
# NMD — Main Figure (FILTERED & REORGANIZED)
#
# v10:
#   * Main B = PTC-to-start distance bins; in-frame AUG and Kozak strength
#     (old C, D) moved to the Supplemental Figure (A, B)
#
# v9:
#   * B, D (and the supplement) use a beeswarm + median/IQR style, distinct
#     from the violin + box style of Figure 3 A-C
#   * Panel E = SHAP cutoff figure (SHAP_cutoff_relativePTClocation.png),
#     spanning both rows on the right; set USE_SHAP_IMAGE <- FALSE to redraw
#     it from the shap_cutoff_*.csv files instead
#
# v8 (reviewer response):
#   * All categorical panels share one style: violin + box + median
#   * Old B (PTC-to-start bins) and old C (in-frame AUG) moved to a separate
#     Supplemental Figure (A, B), written by this same script
#   * Main figure re-lettered: A, B (Kozak), C (rel. PTC scatter),
#     D (PTC tertile), E (SHAP cutoff)
#
# Filter applied throughout:
#   ALLELE.RAT >= 0.35  (same as Figure 3)
#
# Row 1 cohort:
#   A, B = ≥3 exons + exon ≤500 bp                     (df_startprox)
#   C, D = first 200 bp from start (df_first200 / df_aug)
#          NO ≥3-exons or exon-length constraint — only ALLELE.RAT >= 0.35
#   A = Scatter + loess:           PTC.2.start vs NMD (variant-level)
#   B = Raincloud:                 NMD by PTC.2.start bins (variant-level)
#   C = Violin + beeswarm + box:   NMD by in-frame AUG status
#   D = Beeswarm + box:            Kozak strength (with-AUG subset)
#
# Row 2 (Relative PTC location, variant-level):
#   E = Scatter + GAM:             relativePTClocation vs NMD
#   F = Violin + bar + CI:         NMD by PTC tertile
#   G = SHAP cutoff analysis:      Native ggplot reconstruction from CSVs
#                                  (scatter + LOESS, derivative, histogram)
#
# Required files in working directory:
#   TOPMed_stopgain_September25_corrected_readyformodel.csv    (data for Panels A-F)
#   shap_cutoff_curve_data.csv                        (Panel G — variant scatter)
#   shap_cutoff_loess.csv                             (Panel G — LOESS smooth)
#   shap_cutoff_markers.csv                           (Panel G — markers)
#
# Generate the SHAP CSVs once with:
#   python prepare_shap_cutoff_data.py
#
# Removed from this figure (moved to NMD_Supplemental_Figure3.R):
#   - CDS length × PTC location stratification (was G, H)
#   - Exon count × PTC location stratification (was I, J)
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
# library(ggdist)    # v8: not needed (no raincloud panels any more)
# library(gghalves)   # no longer needed (see left-side point layer)
library(ggbeeswarm)   # v9: beeswarm points in B, D and the supplement
library(ggforce)
library(patchwork)
library(scales)

# Resolve namespace conflicts
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
distinct  <- dplyr::distinct
left_join <- dplyr::left_join
arrange   <- dplyr::arrange

set.seed(42)
# (n_iter removed — no bootstrap in this figure version)

# ══════════════════════════════════════════════════════════════════════════════
# COLOURS & THEME
# ══════════════════════════════════════════════════════════════════════════════
TITLE_COL <- "#2E6DA4"
FONT      <- 36   # was 44 — dropped slightly so axis text fits without overlap

# PTC-to-start distance bins (panels A–D)
pal_bins <- c(
  "0-100"   = "#4CAF7D",
  "100-200" = "#E88C3A",
  "200-300" = "#4A72B0",
  "300-700" = "#3EC9A7",
  ">700"    = "#D05050"
)

# In-frame AUG (panel C)
pal_aug <- c("Without AUG" = "#7B68EE",
             "With AUG"    = "#D4507A")

# Kozak strength (panel D)
pal_kozak <- c("Weak"     = "#3EC9A7",
               "Moderate" = "#4A72B0",
               "Strong"   = "#D05050")

# PTC tertile (panels E–I)
pal_ptc       <- c("Low"    = "#5B8A3C",
                   "Medium" = "#4A72B0",
                   "High"   = "#A0403A")
pal_ptc_light <- c("Low"    = "#A8CC88",
                   "Medium" = "#92B4E0",
                   "High"   = "#D89090")

# CDS length (panel H)
pal_cds <- c("Short"  = "#E67E22",
             "Medium" = "#2980B9",
             "Long"   = "#27AE60")

# Exon count (panels I, J)
pal_exon <- c("2-5"   = "#E74C3C",
              "6-10"  = "#E67E22",
              "11-20" = "#2980B9",
              ">20"   = "#27AE60")

# ── Base themes ───────────────────────────────────────────────────────────────
base_theme <- theme_classic(base_size = FONT) +
  theme(
    legend.position    = "none",
    axis.text.x        = element_text(size = FONT + 2, face = "bold",
                                      color = "grey10", lineheight = 1.1),
    axis.text.y        = element_text(size = FONT + 1, face = "bold",
                                      color = "grey15"),
    axis.title.x       = element_text(size = FONT + 3, face = "bold",
                                      color = "grey10",
                                      margin = margin(t = 12)),
    axis.title.y       = element_text(size = FONT + 3, face = "bold",
                                      color = "grey10",
                                      margin = margin(r = 12)),
    axis.line          = element_line(color = "grey20", linewidth = 1.4),
    axis.ticks         = element_line(color = "grey35", linewidth = 1.0),
    panel.grid.major.y = element_line(color = "grey91", linewidth = 0.65),
    panel.grid.major.x = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    strip.text         = element_text(size = FONT + 1, face = "bold",
                                      color = "grey10",
                                      margin = margin(7, 5, 7, 5)),
    strip.background   = element_rect(fill = "grey94", color = "grey60",
                                      linewidth = 1.1),
    plot.tag           = element_text(size = FONT + 18, face = "bold",
                                      color = "grey10"),
    plot.title         = element_text(size = FONT + 4, face = "bold",
                                      color = TITLE_COL, hjust = 0.5,
                                      margin = margin(b = 8)),
    plot.subtitle      = element_text(size = FONT, face = "bold.italic",
                                      color = "grey40", hjust = 0.5,
                                      margin = margin(b = 6)),
    plot.caption       = element_text(size = FONT - 4, face = "italic",
                                      color = "grey50",
                                      hjust = 0.5, margin = margin(t = 8)),
    plot.margin        = margin(16, 22, 16, 22),
    panel.spacing      = unit(1.2, "lines")
  )

scatter_theme <- base_theme +
  theme(
    panel.grid.major   = element_line(color = "grey91", linewidth = 0.65),
    panel.grid.major.x = element_line(color = "grey91", linewidth = 0.65)
  )

# Shared scale objects
ref05      <- geom_hline(yintercept = 0.5, linetype = "dashed",
                         color = "grey38", linewidth = 1.2)
y_nmd_full <- scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                                   limits = c(-0.04, 1.58), expand = c(0, 0))
y_nmd_scat <- scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                                   limits = c(0, 1.06), expand = c(0, 0))
y_dist     <- scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                                   limits = c(-0.04, 1.52), expand = c(0, 0))
y_scat     <- scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                                   limits = c(-0.02, 1.05), expand = c(0, 0))

# ══════════════════════════════════════════════════════════════════════════════
# LOAD & PREP
# ══════════════════════════════════════════════════════════════════════════════
cat("Loading data...\n")
df_raw <- read.csv("TOPMed_stopgain_September25_corrected_readyformodel.csv",
                   stringsAsFactors = FALSE) %>%
  filter(!is.na(ALLELE.RAT), !is.na(PTC.2.start),
         !is.na(TxName), !is.na(cds_length), !is.na(exon_count),
         !is.na(relativePTClocation),
         ALLELE.RAT >= 0.35) %>%      # variant-level quality filter (matches Figure 3)
  mutate(
    # PTC-to-start distance bins
    ptc_bin = factor(PTC.2.start.binning,
                     levels = c("0-100","100-200","200-300","300-700",">700")),

    # In-frame AUG status (panel C)
    aug_status = case_when(
      kozak_strength %in% c("weak","moderate","strong") ~ "With AUG",
      TRUE ~ "Without AUG"
    ),
    aug_status = factor(aug_status, levels = c("Without AUG","With AUG")),

    # Kozak strength (panel D)
    kozak_grp = case_when(
      kozak_strength == "weak"     ~ "Weak",
      kozak_strength == "moderate" ~ "Moderate",
      kozak_strength == "strong"   ~ "Strong",
      TRUE ~ NA_character_
    ),
    kozak_grp = factor(kozak_grp, levels = c("Weak","Moderate","Strong")),

    # PTC tertile (panels E–I)
    ptc_cat = case_when(
      relativePTClocation <= 0.3750 ~ "Low",
      relativePTClocation <= 0.7144 ~ "Medium",
      TRUE                          ~ "High"
    ),
    ptc_cat = factor(ptc_cat, levels = c("Low","Medium","High")),

    # CDS length groups (panels G, H)
    cds_cat = case_when(
      cds_length <  1289 ~ "Short",
      cds_length <= 2181 ~ "Medium",
      TRUE               ~ "Long"
    ),
    cds_cat = factor(cds_cat, levels = c("Short","Medium","Long")),

    # Exon count groups (panels I, J)
    exon_4grp = case_when(
      exon_count <=  5 ~ "2-5",
      exon_count <= 10 ~ "6-10",
      exon_count <= 20 ~ "11-20",
      TRUE             ~ ">20"
    ),
    exon_4grp = factor(exon_4grp, levels = c("2-5","6-10","11-20",">20"))
  )

# ══════════════════════════════════════════════════════════════════════════════
# START-PROXIMAL COHORT for Panels A–D
#   - transcripts with >= 3 exons (to allow meaningful PTC/start architecture)
#   - PTC sits in an exon <= 500 bp (length.mutated.exon)
#   - This is on top of the global ALLELE.RAT >= 0.35 filter applied above
# ══════════════════════════════════════════════════════════════════════════════
if (!"length.mutated.exon" %in% names(df_raw)) {
  stop("Column 'length.mutated.exon' not found in df_raw.")
}

df_startprox <- df_raw %>%
  filter(exon_count >= 3,
         !is.na(length.mutated.exon),
         length.mutated.exon <= 500)

cat(sprintf("\n[Start-proximal cohort: exon_count >= 3 AND length.mutated.exon <= 500]\n"))
cat(sprintf("Variants in start-proximal cohort: %d (%.1f%% of df_raw)\n",
            nrow(df_startprox),
            100 * nrow(df_startprox) / nrow(df_raw)))

# Subsets for the four Row-1 panels
#   A, B use df_startprox            (≥3 exons + exon ≤500 bp)
#   C, D use df_first200 / df_aug    (df_raw + first 200 bp + optional Kozak)
#                                    NOT built on df_startprox
df_first200 <- df_raw       %>% filter(first.200 == "first 200")
df_aug      <- df_first200 %>% drop_na(kozak_grp)

cat(sprintf("All variants (df_raw):              %d\n", nrow(df_raw)))
cat(sprintf("Start-proximal [Panels A, B]:       %d\n", nrow(df_startprox)))

# Sanity check: what values does the first.200 column take?
if (!"first.200" %in% names(df_raw)) {
  stop("Column 'first.200' not found in df_raw. Available columns:\n  ",
       paste(names(df_raw), collapse = ", "))
}
cat("\n[Diagnostic] Distribution of first.200 in df_raw:\n")
print(table(df_raw$first.200, useNA = "ifany"))

cat(sprintf("first 200 bp cohort [Panels C, D]:  %d  (only ALLELE.RAT >= 0.35 + first.200 == 'first 200')\n",
            nrow(df_first200)))
cat(sprintf("+ AUG annotation [Panel D only]:    %d\n", nrow(df_aug)))

# ══════════════════════════════════════════════════════════════════════════════
# PANEL C COHORT — First 200 bp variants, classified by aug_distance_category.
#                  Variants with aug_distance_category == "no_inframe_AUG" are
#                  the "Without AUG" group; all remaining variants are "With AUG".
#                  Reference target: 59 (No) / 763 (Yes).
# ══════════════════════════════════════════════════════════════════════════════
if (!"aug_distance_category" %in% names(df_raw)) {
  stop("Column 'aug_distance_category' not found in df_raw. Column names available:\n  ",
       paste(names(df_raw), collapse = ", "))
}

cat("\n--- Panel C: aug_distance_category diagnostic ---\n")
cat("Unique values of aug_distance_category in df_raw:\n")
print(table(df_raw$aug_distance_category, useNA = "ifany"))

cat("\nUnique values in the first.200 cohort only:\n")
print(table(df_first200$aug_distance_category, useNA = "ifany"))

# Build panel C dataset.
# Drop rows where aug_distance_category is NA (can't classify).
# Then: "no_inframe_AUG" -> Without AUG ; anything else -> With AUG.
df_panelC <- df_first200 %>%
  filter(!is.na(aug_distance_category)) %>%
  mutate(aug_status = if_else(aug_distance_category == "no_inframe_AUG",
                              "Without AUG", "With AUG"),
         aug_status = factor(aug_status, levels = c("Without AUG","With AUG")))

cat("\n>>> Panel C final counts (target: 59 No / 763 Yes):\n")
print(table(df_panelC$aug_status))
cat(sprintf("    total n = %d\n", nrow(df_panelC)))

# ── Cross-check: alternative reading where 'remove' meant delete entirely ────
# If the instruction meant drop no_inframe_AUG rows and keep only variants with
# an in-frame AUG (single group, no comparison), this is what that cohort size
# would be. Printed purely as a sanity check.
cat(sprintf("\n    (If 'remove' meant delete: with-AUG-only cohort would be n = %d)\n\n",
            df_first200 %>%
              filter(!is.na(aug_distance_category),
                     aug_distance_category != "no_inframe_AUG") %>%
              nrow()))

# ══════════════════════════════════════════════════════════════════════════════
# (Bootstrap block removed — was used only for the old Panels E–J which are
#  either now variant-level (E, F) or have moved to NMD_Supplemental_Figure3.R
#  (CDS/exon stratification).)
# ══════════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════════
# PANEL A — Scatter + loess: PTC.2.start vs NMD (variant-level)
#           Uses df_startprox (≥3 exons, exon ≤500 bp, ALLELE.RAT >= 0.35)
#           NOT restricted to first 200 bp — Panel A's x-axis spans 0-1000 nt
#           x-axis limited to 1000 nt (close-up view of the active range)
# ══════════════════════════════════════════════════════════════════════════════
x_max_A <- 1000
df_A    <- df_startprox %>% filter(PTC.2.start <= x_max_A)
n_rm_A  <- sum(df_startprox$PTC.2.start > x_max_A, na.rm = TRUE)
r_A     <- cor(df_A$PTC.2.start, df_A$ALLELE.RAT, use = "complete.obs",
               method = "spearman")   # <<< v6: Spearman (was default Pearson)

pA <- ggplot(df_A, aes(x = PTC.2.start, y = ALLELE.RAT)) +
  geom_point(color = "#A05070", fill = "#A05070",
             size = 3.6, alpha = 0.45, shape = 16) +   # size was 2.0 (change #1)
  geom_smooth(method = "loess", span = 0.50,
              color = "black", linewidth = 2.5,
              fill = "grey72", alpha = 0.25, se = TRUE) +
  ref05 +
  annotate("text", x = x_max_A * 0.62, y = 0.07,
           label = sprintf("\u03c1 = %.3f", r_A),
           size = 14, fontface = "bold.italic", color = "grey20") +
  annotate("text", x = x_max_A * 0.01, y = 0.98,
           label = paste0("n = ", format(nrow(df_A), big.mark = ",")),
           hjust = 0, size = 14, fontface = "bold", color = "grey20") +
  scale_x_continuous(breaks = seq(0, x_max_A, by = 200),
                     limits = c(0, x_max_A * 1.02),
                     labels = comma, expand = c(0.01, 0)) +
  y_nmd_scat +
  labs(x        = "PTC-to-start distance (nt)",
       y        = "NMD efficiency",
       title    = "PTC-to-Start Distance vs NMD Efficiency",
       tag      = "A") +
  scatter_theme

# ══════════════════════════════════════════════════════════════════════════════
# SHARED CATEGORICAL STYLE (v9) — used for every categorical panel:
#   main Figure 4 B (Kozak) and D (PTC tertile),
#   Supplemental Figure A (PTC-to-start bins) and B (in-frame AUG)
#   beeswarm of all variants + median (large dot) and IQR (thick bar);
#   deliberately different from Figure 3 A–C (violin + box),
#   same y-axis, "n = …" labels, Wilcoxon brackets (pairs with >=3 obs only)
# ══════════════════════════════════════════════════════════════════════════════
make_violin_panel <- function(d, xvar, pal, pairs, bracket_y,
                              xlab, title, tag, x_angle = 0, n_size = 11) {
  d <- d %>% drop_na(all_of(xvar))

  n_lab <- d %>% count(across(all_of(xvar))) %>%
    mutate(label = paste0("n = ", format(n, big.mark = ",", trim = TRUE)))

  summ <- d %>% group_by(across(all_of(xvar))) %>%
    summarise(med = median(ALLELE.RAT),
              q25 = quantile(ALLELE.RAT, 0.25),
              q75 = quantile(ALLELE.RAT, 0.75), .groups = "drop")

  viable <- as.character(n_lab[[xvar]][n_lab$n >= 3])
  valid  <- Filter(function(pr) all(pr %in% viable), pairs)

  st <- NULL
  if (length(valid) > 0) {
    st <- tryCatch({
      d %>%
        wilcox_test(as.formula(paste("ALLELE.RAT ~", xvar)),
                    comparisons = valid) %>%
        add_significance("p") %>%
        add_xy_position(x = xvar) %>%
        mutate(y.position = rep_len(bracket_y, n()),
               !!xvar := group1)
    }, error = function(e) {
      cat(sprintf("[%s] Wilcoxon failed (%s); skipping brackets.\n",
                  tag, conditionMessage(e)))
      NULL
    })
  }

  p <- ggplot(d, aes(x = .data[[xvar]], y = ALLELE.RAT)) +
    # every variant as a beeswarm point ...
    geom_quasirandom(aes(color = .data[[xvar]]), size = 2.4, alpha = 0.35,
                     shape = 16, width = 0.36) +
    # ... with the median (large dot) and interquartile range (thick bar)
    geom_linerange(data = summ, aes(x = .data[[xvar]], ymin = q25, ymax = q75),
                   inherit.aes = FALSE, linewidth = 4.2, color = "grey10") +
    geom_point(data = summ, aes(x = .data[[xvar]], y = med, fill = .data[[xvar]]),
               inherit.aes = FALSE, shape = 21, size = 10,
               color = "grey10", stroke = 2.2) +
    ref05 +
    geom_text(data = n_lab, aes(x = .data[[xvar]], y = 1.48, label = label),
              inherit.aes = FALSE, size = n_size,
              fontface = "bold", color = "grey15") +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) +
    y_nmd_full +
    labs(x = xlab, y = "NMD efficiency", title = title, tag = tag) +
    base_theme +
    theme(axis.text.x = element_text(size = FONT - 2, face = "bold",
                                     color = "grey10", angle = x_angle,
                                     hjust = if (x_angle > 0) 1 else 0.5,
                                     vjust = if (x_angle > 0) 1 else 0.5))

  if (!is.null(st) && nrow(st) > 0)
    p <- p + stat_pvalue_manual(st, label = "p.signif",
                                tip.length = 0.012, bracket.size = 1.0,
                                size = 14, color = "grey20")
  p
}

# ── Main B — NMD by PTC-to-start distance bins ───────────────────────────────
#     df_startprox (>=3 exons, exon <=500 bp, ALLELE.RAT >= 0.35)
df_B <- df_startprox %>% drop_na(ptc_bin)
cat("\n[Main B] PTC.2.start bin counts in df_startprox:\n")
print(df_B %>% count(ptc_bin))

pB <- make_violin_panel(
  df_B, "ptc_bin", pal_bins,
  pairs = list(c("0-100","100-200"), c("100-200","200-300"),
               c("200-300","300-700"), c("300-700",">700")),
  bracket_y = c(1.08, 1.18, 1.28, 1.18),
  xlab = "PTC-to-start distance (nt)",
  title = "NMD by PTC-to-Start Distance Bins", tag = "B",
  x_angle = 30, n_size = 9.5)          # 5 groups: tilt labels, smaller n

# ── Supplemental A — NMD by in-frame AUG (was main Panel C) ──────────────────
pS_A <- make_violin_panel(
  df_panelC, "aug_status", pal_aug,
  pairs = list(c("Without AUG", "With AUG")), bracket_y = 1.20,
  xlab = "In-frame AUG status",
  title = "NMD Efficiency by In-Frame AUG", tag = "A")

# ── Supplemental B — Kozak strength (was main Panel D) ───────────────────────
pS_B <- make_violin_panel(
  df_aug, "kozak_grp", pal_kozak,
  pairs = list(c("Weak","Moderate"), c("Moderate","Strong"), c("Weak","Strong")),
  bracket_y = c(1.10, 1.22, 1.34),
  xlab = "Kozak strength", title = "Kozak Strength", tag = "B")

# ══════════════════════════════════════════════════════════════════════════════
# PANEL E — Scatter + GAM: relativePTClocation vs NMD (VARIANT-LEVEL)
#           One point per variant in df_raw (filtered to ALLELE.RAT >= 0.35).
# ══════════════════════════════════════════════════════════════════════════════
df_E_var <- df_raw %>% drop_na(relativePTClocation, ALLELE.RAT)
r_E_var  <- cor(df_E_var$relativePTClocation, df_E_var$ALLELE.RAT,
                use = "complete.obs",
                method = "spearman")   # <<< v6: Spearman (was default Pearson)

pE <- ggplot(df_E_var,
             aes(x = relativePTClocation, y = ALLELE.RAT)) +
  geom_point(color = "#4DA899", fill = "#4DA899",
             size = 3.0, alpha = 0.35, shape = 21, stroke = 0.25) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              color = "black", linewidth = 2.4,
              fill = "grey70", alpha = 0.22, se = TRUE) +
  ref05 +
  annotate("text", x = 0.04, y = 0.96,
           label = paste0("n = ", format(nrow(df_E_var), big.mark = ",")),
           hjust = 0, size = 14, fontface = "bold", color = "grey20") +
  annotate("text", x = 0.62, y = 0.07,
           label = sprintf("\u03c1 = %.3f", r_E_var),
           size = 14, fontface = "bold.italic", color = "grey20") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                     limits = c(0, 1), expand = c(0.01, 0)) +
  y_scat +
  labs(x       = "Relative PTC location",
       y       = "NMD efficiency",
       title   = "PTC Location vs NMD Efficiency",
       tag     = "C") +
  scatter_theme

# ══════════════════════════════════════════════════════════════════════════════
# Main D — NMD by relative PTC location tertile (was Panel F; now shared style)
# ══════════════════════════════════════════════════════════════════════════════
df_F_var <- df_raw %>% drop_na(ptc_cat, ALLELE.RAT)

pD <- make_violin_panel(
  df_F_var, "ptc_cat", pal_ptc,
  pairs = list(c("Low","Medium"), c("Medium","High"), c("Low","High")),
  bracket_y = c(1.10, 1.20, 1.30),
  xlab = "Relative PTC location (tertile)",
  title = "NMD by Relative PTC Location", tag = "D")

# ══════════════════════════════════════════════════════════════════════════════
# PANEL G — SHAP cutoff plot for relativePTClocation (native ggplot)
#
# Reads three CSVs produced by prepare_shap_cutoff_data.py:
#   shap_cutoff_curve_data.csv  — variant-level scatter (relativePTClocation, shap_value)
#   shap_cutoff_loess.csv       — LOESS smooth + derivative (x_grid, y_smooth, dy_dx)
#   shap_cutoff_markers.csv     — markers (zero crossings + rate-of-change peaks)
#
# To regenerate the CSVs:
#   python prepare_shap_cutoff_data.py
# ══════════════════════════════════════════════════════════════════════════════
shap_files <- c("shap_cutoff_curve_data.csv",
                "shap_cutoff_loess.csv",
                "shap_cutoff_markers.csv")

# Panel E source:
#   USE_SHAP_IMAGE = TRUE  -> place the finished SHAP figure (PNG) as Panel E
#   USE_SHAP_IMAGE = FALSE -> redraw it in ggplot from the shap_cutoff_*.csv files
USE_SHAP_IMAGE <- FALSE   # TRUE = use the PNG instead
G_IMAGE <- "SHAP_cutoff_relativePTClocation.png"

if ((USE_SHAP_IMAGE || !all(file.exists(shap_files))) && file.exists(G_IMAGE)) {
  g_img <- png::readPNG(G_IMAGE)
  pG <- wrap_elements(full = grid::rasterGrob(g_img, interpolate = TRUE)) +
    labs(tag = "E") +
    theme(plot.tag = element_text(size = FONT + 18, face = "bold",
                                  color = "grey10"))
} else if (any(!file.exists(shap_files))) {
  warning("Missing SHAP cutoff CSVs (run `python prepare_shap_cutoff_data.py` first); ",
          "Panel G will render as a placeholder.")
  pG <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0("Panel E placeholder\n",
                            "Run prepare_shap_cutoff_data.py to\n",
                            "generate shap_cutoff_*.csv"),
             size = 13, color = "grey30") +
    theme_void() +
    theme(plot.tag = element_text(size = FONT + 14, face = "bold")) +
    labs(tag = "E")
} else {
  shap_scatter <- read.csv("shap_cutoff_curve_data.csv", stringsAsFactors = FALSE)
  shap_loess   <- read.csv("shap_cutoff_loess.csv",      stringsAsFactors = FALSE)
  shap_marks   <- read.csv("shap_cutoff_markers.csv",    stringsAsFactors = FALSE)

  zero_marks <- shap_marks %>% filter(marker_type == "zero_crossing")
  peak_marks <- shap_marks %>% filter(marker_type == "rate_peak")

  # Color helpers
  zero_col <- function(direction) {
    ifelse(direction == "pos_to_neg", "#3b78b5",   # blue
                                       "#3aa86b")  # green
  }
  peak_col <- function(direction) {
    ifelse(direction == "decreasing", "#3b78b5",   # blue (descending side)
                                       "#cc3333")  # red (ascending side)
  }

  # ── Top sub-panel: scatter + LOESS + zero-crossing markers ───────────────
  pG_top <- ggplot() +
    geom_point(data = shap_scatter,
               aes(x = relativePTClocation, y = shap_value),
               color = "grey55", alpha = 0.25, size = 2.2, shape = 16) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.7) +
    geom_line(data = shap_loess,
              aes(x = x_grid, y = y_smooth),
              color = "#cc3333", linewidth = 3.4) +

    # vertical lines at zero crossings
    geom_segment(data = zero_marks,
                 aes(x = x, xend = x, y = -Inf, yend = Inf,
                     color = label),
                 linewidth = 2.0, show.legend = FALSE) +

    # labeled markers near top
    geom_label(data = zero_marks,
               aes(x = x, y = max(shap_scatter$shap_value) + 0.10,
                   label = sprintf("%.2f", x),
                   color = label),
               size = 15, fontface = "bold",
               fill = alpha("white", 0.92),
               label.size = 0.9,
               show.legend = FALSE) +

    scale_color_manual(values = c("pos_to_neg" = "#3b78b5",
                                   "neg_to_pos" = "#3aa86b")) +

    scale_x_continuous(breaks = seq(0, 1, 0.2),
                       limits = c(0, 1), expand = c(0.005, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
    labs(x = NULL, y = "SHAP value",
         title = "SHAP cutoff analysis: Relative PTC Location",
         tag = "E") +
    base_theme +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = FONT + 4, face = "bold",
                                       margin = margin(r = 6)),
          axis.text.y  = element_text(size = FONT + 2),
          plot.title   = element_text(size = FONT + 6, face = "bold",
                                       color = TITLE_COL, hjust = 0.5,
                                       margin = margin(b = 8)))

  # ── Middle sub-panel: derivative d(SHAP)/d(feature) with rate-of-change peaks
  # Pull dy values at each peak for plotting markers
  peak_marks <- peak_marks %>%
    mutate(dy_val = approx(shap_loess$x_grid, shap_loess$dy_dx, xout = x)$y)

  pG_mid <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.7) +
    geom_line(data = shap_loess,
              aes(x = x_grid, y = dy_dx),
              color = "#9466bb", linewidth = 2.4) +

    # vertical dotted lines at peaks
    geom_segment(data = peak_marks,
                 aes(x = x, xend = x, y = -Inf, yend = dy_val,
                     color = label),
                 linetype = "dotted", linewidth = 1.4,
                 show.legend = FALSE) +
    geom_point(data = peak_marks,
               aes(x = x, y = dy_val, fill = label),
               shape = 21, size = 7.0, stroke = 1.5,
               color = "grey15", show.legend = FALSE) +
    geom_text(data = peak_marks,
              aes(x = x + 0.035 * ifelse(label == "increasing", -1, 1),
                  y = dy_val * 0.85,
                  hjust = ifelse(label == "increasing", 1, 0),
                  label = sprintf("%.2f", x),
                  color = label),
              size = 15, fontface = "bold", show.legend = FALSE) +

    scale_color_manual(values = c("decreasing" = "#3b78b5",
                                   "increasing" = "#cc3333")) +
    scale_fill_manual( values = c("decreasing" = "#3b78b5",
                                   "increasing" = "#cc3333")) +

    scale_x_continuous(breaks = seq(0, 1, 0.2),
                       limits = c(0, 1), expand = c(0.005, 0)) +
    scale_y_continuous(expand = expansion(mult = 0.22)) +
    labs(x = NULL, y = "d(SHAP)/dx") +
    base_theme +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = FONT + 4, face = "bold",
                                       margin = margin(r = 6)),
          axis.text.y  = element_text(size = FONT + 2))

  # ── Bottom sub-panel: variant histogram with same vertical markers ─────────
  pG_bot <- ggplot(shap_scatter, aes(x = relativePTClocation)) +
    geom_histogram(bins = 50, fill = "grey75",
                   color = "grey25", linewidth = 0.5) +

    geom_segment(data = zero_marks,
                 aes(x = x, xend = x, y = -Inf, yend = Inf,
                     color = label),
                 linewidth = 1.6, show.legend = FALSE) +

    scale_color_manual(values = c("pos_to_neg" = "#3b78b5",
                                   "neg_to_pos" = "#3aa86b")) +

    scale_x_continuous(breaks = seq(0, 1, 0.2),
                       limits = c(0, 1), expand = c(0.005, 0)) +
    labs(x = "Relative PTC location", y = "Variants") +
    base_theme +
    theme(axis.text.x  = element_text(size = FONT + 4, face = "bold"),
          axis.title.x = element_text(size = FONT + 6, face = "bold")) +
    theme(axis.title.y = element_text(size = FONT + 4, face = "bold",
                                       margin = margin(r = 6)),
          axis.text.y  = element_text(size = FONT + 2))

  # Stack the three sub-panels (each keeps its own y-axis; patchwork
  # auto-aligns the plot panel positions vertically)
  pG <- pG_top / pG_mid / pG_bot +
    plot_layout(heights = c(2.2, 1.2, 0.9))
}

# ══════════════════════════════════════════════════════════════════════════════
# COMBINE — MAIN FIGURE 4 (v8)
#   Row 1:  A (PTC-to-start scatter) | B (PTC-to-start bins) | E (SHAP cutoff,
#   Row 2:  C (relative PTC scatter) | D (PTC tertile)       |   spans both rows)
#
#   Panel mapping from v7:  A→A, B→B, E→C, F→D, G→E
#   Old C (in-frame AUG) and old D (Kozak strength) → Supplemental Figure
# ══════════════════════════════════════════════════════════════════════════════
# Panel E (SHAP, three stacked sub-panels) spans both rows on the right so
# each sub-panel is tall enough for the larger fonts.
fig4_design <- "
AAAABBBBEEEEE
CCCCDDDDEEEEE
"
combined <- wrap_plots(A = pA, B = pB, C = pE, D = pD, E = pG,
                       design = fig4_design) &
  theme(plot.margin = margin(12, 40, 12, 16))   # extra right margin: last x tick not clipped

ggsave("Figure4_revised.pdf", combined, device = cairo_pdf,
       width = 44, height = 24, dpi = 300, limitsize = FALSE)
ggsave("Figure4_revised.png", combined,
       width = 44, height = 24, dpi = 300, limitsize = FALSE)

print(combined)
cat("Saved: Figure4_revised.pdf / .png\n")

# ══════════════════════════════════════════════════════════════════════════════
# SUPPLEMENTAL FIGURE — start-proximal PTCs (same violin style as Figure 4)
#   A: NMD by in-frame AUG status          (first 200 bp cohort)
#   B: NMD by Kozak strength               (first 200 bp + in-frame AUG)
# ══════════════════════════════════════════════════════════════════════════════
supp <- (pS_A | pS_B) +
  plot_layout(widths = c(1, 1)) &
  theme(plot.margin = margin(12, 16, 12, 16))

ggsave("Supplemental_Figure_StartProximal.pdf", supp, device = cairo_pdf,
       width = 24, height = 12, dpi = 300, limitsize = FALSE)
ggsave("Supplemental_Figure_StartProximal.png", supp,
       width = 24, height = 12, dpi = 300, limitsize = FALSE)
cat("Saved: Supplemental_Figure_StartProximal.pdf / .png\n")
