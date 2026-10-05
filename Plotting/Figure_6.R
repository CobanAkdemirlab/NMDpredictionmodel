# ══════════════════════════════════════════════════════════════════════════════
# NMD — PTC-to-EJC Region Features vs NMD Efficiency  (v11, variant-level)
#
# CHANGES vs v10 (Reviewer 1, specific comment 8: labels in panel B cut off /
# overlapping between neighbouring facets):
#   - Panel B facets: wider spacing between facets (2.2 -> 5 lines), extra
#     right-hand padding, and the last x tick pulled inside each facet
#     (100->80, 200->150, 500->400, 1000->800) so end labels of adjacent
#     facets no longer collide
#   - Panel B statistics: rho and p on separate lines, slightly smaller, so
#     they fit inside each facet instead of running past its edge
#   - The standalone codon-optimality SHAP plot (not part of this figure) is
#     now OFF by default (BUILD_SHAP_STANDALONE <- FALSE), so the figure no
#     longer requires Python or the SHAP pickle
#
# CHANGES vs v7:
#   - PTC-variant filter:  ALLELE.RAT >= 0.35 (matches Figures 3 and 4)
#   - Panel D: dropped the "Whole CDS" facet — keeps only PTC ± 100 nt
#   - Larger base font, thicker lines, bigger panel-letter tags
#
# Inherited from v7:
#   - Variant-level analysis (NO bootstrap; ALLELE.RAT used directly)
#   - GLOBAL subset applied to every panel:
#         last.EJC == "upstream"   AND
#         first.200 != "first 200"   (i.e. PTC outside first 200 bp of CDS)
#   - Panel A: NMD efficiency by PTC-containing exon length bins
#         Short <100 | Medium 100–200 | Long 200–500 | Very long >500
#
#   A = Hex+LOESS      NMD efficiency vs PTC-to-EJC distance (overall)
#   B = Hex+LOESS      Same, faceted by PTC-exon length (4 bins)
#   C = Bar            PhastCons sequential p-values
#   D = Lollipop       Top10 inhibiting + top10 promoting RBP motifs
#                       (PTC-to-EJC region)
#   E = Lollipop       Bonferroni-significant RBP motifs — PTC±100 nt
#   F = Lollipop       Bonferroni-significant RBP motifs — EJC±100 nt
#
#   (The SHAP codon-optimality cutoff scatter — Low <0.4 / Medium 0.4-0.6 /
#    High >0.6, PTC ± 100 nt — is built below but is NOT part of this figure;
#    it is saved standalone and omitted from the merged layout. Panels are
#    lettered A–F with no gap.)
#
# Required file: TOPMed_stopgain_September25_corrected_readyformodel.csv
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggdist)
library(gghalves)
library(ggbeeswarm)
library(ggforce)
library(patchwork)
library(scales)
# cowplot / magick are only needed for the optional standalone SHAP plot
if (requireNamespace("cowplot", quietly = TRUE)) library(cowplot)
if (requireNamespace("magick",  quietly = TRUE)) library(magick)

select    <- dplyr::select;   filter    <- dplyr::filter
mutate    <- dplyr::mutate;   summarise <- dplyr::summarise
group_by  <- dplyr::group_by; distinct  <- dplyr::distinct
left_join <- dplyr::left_join; arrange  <- dplyr::arrange

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

# ── Colours & theme ───────────────────────────────────────────────────────────
TITLE_COL <- "#2E6DA4"
FONT      <- 96        # bumped from 80
RBP_MOTIF_SIZE <- 64   # bumped from 56

base_theme <- theme_classic(base_size = FONT) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(size = FONT,     face = "bold",
                                    color = "grey5",  lineheight = 1.1),
    axis.text.y      = element_text(size = FONT,     face = "bold",
                                    color = "grey5"),
    axis.title.x     = element_text(size = FONT + 8, face = "bold",
                                    color = "grey5",  margin = margin(t = 26)),
    axis.title.y     = element_text(size = FONT + 8, face = "bold",
                                    color = "grey5",  margin = margin(r = 26)),
    axis.line        = element_line(color = "grey10", linewidth = 4.4),  # was 3.4
    axis.ticks       = element_line(color = "grey15", linewidth = 3.4),  # was 2.6
    axis.ticks.length = unit(0.72, "cm"),
    panel.grid.major = element_line(color = "grey88", linewidth = 1.8),  # was 1.4
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    strip.text       = element_text(size = FONT + 4, face = "bold",
                                    color = "grey5",  margin = margin(16, 6, 16, 6)),
    strip.background = element_rect(fill = "grey92", color = "grey40",
                                    linewidth = 2.4),
    plot.tag         = element_text(size = FONT + 50, face = "bold",      # was +40
                                    color = "grey5"),
    plot.title       = element_text(size = FONT + 12, face = "bold",      # was +10
                                    color = TITLE_COL, hjust = 0.5,
                                    margin = margin(b = 24)),
    plot.margin      = margin(28, 32, 28, 32),
    panel.spacing    = unit(2.4, "lines")
  )

ref05      <- geom_hline(yintercept = 0.5, linetype = "dashed",
                         color = "grey30", linewidth = 3.4)
y_nmd_full <- scale_y_continuous(breaks = c(0, .25, .5, .75, 1),
                                   limits = c(-.04, 1.56), expand = c(0, 0))
y_nmd_bar  <- scale_y_continuous(breaks = seq(0, 1, .25),
                                   limits = c(0, 1.05),   expand = c(0, 0))

seq_comp <- function(lvls)
  lapply(seq_len(length(lvls) - 1), function(i) c(lvls[i], lvls[i + 1]))

# ══════════════════════════════════════════════════════════════════════════════
# LOAD & GLOBAL SUBSET
# ══════════════════════════════════════════════════════════════════════════════
df_raw <- read.csv("TOPMed_stopgain_September25_corrected_readyformodel.csv",
                   stringsAsFactors = FALSE) %>%
  filter(!is.na(ALLELE.RAT), ALLELE.RAT >= 0.35, !is.na(TxName))

cat(sprintf("Variants after ALLELE.RAT >= 0.35 filter: %s\n",
            format(nrow(df_raw), big.mark = ",")))

# Identify RBP motif column sets
rbp_cols_ejc    <- grep("^ptc_to_ejc\\.",     names(df_raw), value = TRUE)
rbp_cols_ptc100 <- grep("^ptcpm100\\.",        names(df_raw), value = TRUE)
rbp_cols_ejc100 <- grep("^ejcpm100\\.",        names(df_raw), value = TRUE)
if (length(rbp_cols_ptc100) == 0)
  rbp_cols_ptc100 <- grep("^ptc_pm100\\.|^ptc\\.pm100\\.", names(df_raw), value = TRUE)
if (length(rbp_cols_ejc100) == 0)
  rbp_cols_ejc100 <- grep("^ejc_pm100\\.|^ejc\\.pm100\\.", names(df_raw), value = TRUE)

cat(sprintf("RBP columns — PTC-to-EJC: %d  |  PTC\u00b1100: %d  |  EJC\u00b1100: %d\n",
            length(rbp_cols_ejc), length(rbp_cols_ptc100), length(rbp_cols_ejc100)))

# Categorical features used by panels B and C
df_raw <- df_raw %>%
  mutate(
    phast_cat = case_when(
      phastcons_ptc_to_ejc_median == 0   ~ "None\n(= 0)",
      phastcons_ptc_to_ejc_median <  0.5 ~ "Low\n(0\u20130.5)",
      phastcons_ptc_to_ejc_median <  1.0 ~ "High\n(0.5\u20131)",
      phastcons_ptc_to_ejc_median == 1.0 ~ "Max\n(= 1)"
    ),
    phast_cat = factor(phast_cat,
                       levels = c("None\n(= 0)", "Low\n(0\u20130.5)",
                                  "High\n(0.5\u20131)", "Max\n(= 1)"))
  )

# ── GLOBAL VARIANT-LEVEL SUBSET ──────────────────────────────────────────────
# Upstream PTC (last.EJC == "upstream") AND not in first 200 bp of CDS
df_var <- df_raw %>%
  filter(last.EJC == "upstream",
         first.200 != "first 200")

cat(sprintf("\nGlobal subset: %s of %s variants kept (upstream PTC, outside first 200 bp)\n",
            format(nrow(df_var), big.mark = ","),
            format(nrow(df_raw), big.mark = ",")))

# ══════════════════════════════════════════════════════════════════════════════
# PANEL A — NMD efficiency vs PTC-to-EJC distance (continuous)
# ══════════════════════════════════════════════════════════════════════════════
dfAB <- df_var %>%
  filter(!is.na(ALLELE.RAT), !is.na(PTC.2.EJC),
         !is.na(length.mutated.exon), PTC.2.EJC > 0) %>%
  mutate(exon_grp = cut(length.mutated.exon,
                         breaks = c(-Inf, 100, 200, 500, Inf),
                         labels = c("Short (<100)",
                                    "Medium (100\u2013200)",
                                    "Long (200\u2013500)",
                                    "Very long (>500)"),
                         right  = FALSE)) %>%
  drop_na(exon_grp)

# Cap PTC.2.EJC at 99th percentile for display only (LOESS uses full data range
# but axis is windowed so a few outliers don't compress the plot)
p99_AB <- quantile(dfAB$PTC.2.EJC, 0.99, na.rm = TRUE)

# Panel A — overall
sp_A <- cor.test(dfAB$PTC.2.EJC, dfAB$ALLELE.RAT, method = "spearman",
                 exact = FALSE)
spearman_label_A <- sprintf("Spearman \u03c1 = %.3f   p = %s",
                            sp_A$estimate,
                            format.pval(sp_A$p.value, digits = 2, eps = 1e-300))

# Plotmath version of the same label (renders rho via R's symbol engine, so it
# does not depend on the text font carrying a Greek glyph). parse = TRUE below.
spearman_label_A_plot <- sprintf('paste("Spearman ", rho, " = %.3f,    ", italic(p), " = %s")',
                                 sp_A$estimate,
                                 format.pval(sp_A$p.value, digits = 2, eps = 1e-300))

cat(sprintf("Panel A: %s variants  |  %s\n",
            format(nrow(dfAB), big.mark = ","), spearman_label_A))

pA <- ggplot(dfAB, aes(x = PTC.2.EJC, y = ALLELE.RAT)) +
  geom_hex(bins = 42, color = "grey85", linewidth = 0.4) +
  geom_smooth(method = "loess", se = TRUE, span = 0.45,
              color = "#C0392B", fill = "#C0392B",
              linewidth = 7.2, alpha = 0.24) +
  ref05 +
  scale_fill_gradientn(colors = c("#F4F8FB", "#9CC2E5", "#2E6DA4", "#10314D"),
                        trans  = "log10",
                        name   = "Variants") +
  scale_x_continuous(limits = c(0, p99_AB),
                     expand = expansion(mult = c(0.01, 0.03)),
                     labels = comma) +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1),
                     limits = c(-.04, 1.04), expand = c(0, 0)) +
  annotate("text", x = p99_AB * 0.97, y = 0.06,
           label = spearman_label_A_plot, parse = TRUE,
           hjust = 1, vjust = 0, size = FONT/3.6,
           fontface = "bold", color = "grey10") +
  annotate("text", x = p99_AB * 0.97, y = 1.00,
           label = sprintf("n = %s", format(nrow(dfAB), big.mark = ",")),
           hjust = 1, vjust = 1, size = FONT/3.6,
           fontface = "bold", color = "grey10") +
  labs(x       = "PTC-to-EJC distance (nt)",
       y       = "NMD efficiency",
       title   = "PTC-to-EJC Distance",
       tag     = "A") +
  base_theme +
  theme(legend.position = "none")

# ══════════════════════════════════════════════════════════════════════════════
# PANEL B — same relationship, faceted by PTC-exon length
#   Per-facet x-axis caps:
#     Short (<100)         -> 0–100
#     Medium (100–200)     -> 0–200
#     Long (200–500)       -> 0–500
#     Very long (>500)     -> 0–1,000
# ══════════════════════════════════════════════════════════════════════════════
EXON_FILL <- c("Short (<100)"          = "#4DA3D9",
               "Medium (100\u2013200)" = "#F9C846",
               "Long (200\u2013500)"   = "#F08A3E",
               "Very long (>500)"      = "#C0392B")

# Per-facet x upper bounds (in nt)
EXON_XMAX <- c("Short (<100)"          =   100,
               "Medium (100\u2013200)" =   200,
               "Long (200\u2013500)"   =   500,
               "Very long (>500)"      = 1000)

# Restrict each facet's data to its visible window so geom_smooth doesn't
# fit on points beyond the displayed range
dfB <- dfAB %>%
  mutate(.xmax = EXON_XMAX[as.character(exon_grp)]) %>%
  filter(PTC.2.EJC <= .xmax)

# "Anchor" frame: invisible points at the desired (x, y) extremes per facet
# so scales = "free_x" honours per-facet caps consistently
anchor_B <- tibble(
  exon_grp   = factor(names(EXON_XMAX), levels = names(EXON_XMAX)),
  PTC.2.EJC  = EXON_XMAX,
  ALLELE.RAT = NA_real_
)

# Per-facet Spearman + n  (computed on the windowed data, matching what is shown)
facet_stats <- dfB %>%
  group_by(exon_grp) %>%
  summarise(rho   = cor(PTC.2.EJC, ALLELE.RAT, method = "spearman",
                         use = "complete.obs"),
            p     = cor.test(PTC.2.EJC, ALLELE.RAT, method = "spearman",
                              exact = FALSE)$p.value,
            n     = n(),
            .groups = "drop") %>%
  mutate(xmax      = EXON_XMAX[as.character(exon_grp)],
         label_rho = sprintf('paste(rho, " = %.3f")', rho),
         label_p   = sprintf('paste(italic(p), " = %s")',
                              format.pval(p, digits = 2, eps = 1e-300)),
         label_n   = sprintf("n = %s", prettyNum(n, big.mark = ",")))

# Choose tick breaks per-facet so they look balanced
# Choose tick breaks per-facet so labels don't collide. Each facet is narrow,
# so we use 3-4 well-spaced breaks instead of fine-grained ones.
B_breaks <- function(limits) {
  ub <- max(limits, na.rm = TRUE)
  if (ub <=  110) c(0, 40, 80)
  else if (ub <=  220) c(0, 75, 150)
  else if (ub <=  550) c(0, 200, 400)
  else                 c(0, 400, 800)
}

pB <- ggplot(dfB, aes(x = PTC.2.EJC, y = ALLELE.RAT, color = exon_grp)) +
  geom_blank(data = anchor_B, inherit.aes = FALSE,
             aes(x = PTC.2.EJC, y = 0.5)) +
  geom_point(alpha = 0.42, size = 9.0, shape = 16) +
  geom_smooth(aes(fill = exon_grp), method = "loess", se = TRUE, span = 0.55,
              linewidth = 7.6, alpha = 0.26) +
  ref05 +
  geom_text(data = facet_stats,
            aes(x = xmax * 0.95, y = 0.17, label = label_rho),
            inherit.aes = FALSE, hjust = 1, vjust = 0, parse = TRUE,
            size = FONT/3.4, fontface = "bold", color = "grey10") +
  geom_text(data = facet_stats,
            aes(x = xmax * 0.95, y = 0.04, label = label_p),
            inherit.aes = FALSE, hjust = 1, vjust = 0, parse = TRUE,
            size = FONT/3.4, fontface = "bold", color = "grey10") +
  geom_text(data = facet_stats,
            aes(x = xmax * 0.95, y = 1.00, label = label_n),
            inherit.aes = FALSE, hjust = 1, vjust = 1,
            size = FONT/3.4, fontface = "bold", color = "grey10") +
  scale_color_manual(values = EXON_FILL) +
  scale_fill_manual(values  = EXON_FILL) +
  scale_x_continuous(breaks = B_breaks,
                     expand = expansion(mult = c(0.02, 0.06)),
                     labels = comma) +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1),
                     limits = c(-.04, 1.04), expand = c(0, 0)) +
  facet_wrap(~ exon_grp, nrow = 1, scales = "free_x") +
  labs(x       = "PTC-to-EJC distance (nt)",
       y       = "NMD efficiency",
       title   = "PTC-to-EJC Distance by PTC-Exon Length",
       tag     = "B") +
  base_theme +
  theme(panel.spacing = unit(5, "lines"))

# ══════════════════════════════════════════════════════════════════════════════
# PANEL C — BAR ± CI: PhastCons sequential p-values
# ══════════════════════════════════════════════════════════════════════════════
dfB   <- df_var %>% drop_na(phast_cat, ALLELE.RAT)
bar_B <- dfB %>%
  group_by(phast_cat) %>%
  summarise(mean_NMD = mean(ALLELE.RAT, na.rm = TRUE),
            se       = sd(ALLELE.RAT, na.rm = TRUE) / sqrt(n()),
            ci_lo    = mean_NMD - 1.96 * se,
            ci_hi    = mean_NMD + 1.96 * se,
            n        = n(), .groups = "drop")

nB <- dfB %>% dplyr::count(phast_cat) %>% mutate(label = format(n, big.mark = ","))

stat_B <- dfB %>%
  wilcox_test(ALLELE.RAT ~ phast_cat,
              comparisons = list(c("None\n(= 0)", "Low\n(0\u20130.5)"),
                                 c("Low\n(0\u20130.5)", "High\n(0.5\u20131)"),
                                 c("High\n(0.5\u20131)", "Max\n(= 1)"))) %>%
  add_significance("p") %>%
  add_xy_position(x = "phast_cat") %>%
  mutate(y.position = c(0.80, 0.88, 0.80), phast_cat = group1)

pal4_B <- c("None\n(= 0)"    = "#BDC3C7",
             "Low\n(0\u20130.5)" = "#85C1E9",
             "High\n(0.5\u20131)"= "#2E86C1",
             "Max\n(= 1)"    = "#1B4F72")

pC <- ggplot(bar_B, aes(x = phast_cat, y = mean_NMD, fill = phast_cat)) +
  geom_col(alpha = 0.85, color = "grey20", linewidth = 2.0, width = 0.72) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.25, linewidth = 3.6, color = "grey15") +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "grey38", linewidth = 3.0) +
  geom_text(aes(label = sprintf("%.3f", mean_NMD), y = ci_hi + 0.015),
            size = 22, fontface = "bold", color = "grey15") +
  geom_text(data = nB, aes(x = phast_cat, y = 0.03, label = label),
            inherit.aes = FALSE, size = 20, fontface = "bold", color = "white") +
  stat_pvalue_manual(stat_B, label = "p.signif", tip.length = 0.016,
                     bracket.size = 2.4, size = 28, color = "grey20") +
  scale_fill_manual(values = pal4_B) +
  y_nmd_bar +
  labs(x       = "PhastCons conservation category",
       y       = "Mean NMD efficiency",
       title   = "PTC-to-EJC Conservation (PhastCons)",
       tag     = "C") +
  base_theme

# ══════════════════════════════════════════════════════════════════════════════
# PANEL D — SHAP CUTOFF SCATTER (codon optimality, PTC ± 100 nt)
#   Embeds the SHAP cutoff-analysis subplot — LOESS over scatter, with the
#   two pooled bootstrap zero-crossings (~+0.41 and ~+0.60) annotated.
#   The image is rebuilt from `shap_values_for_sharing.pkl` on every run
#   (skipped if the PNG is newer than the pickle, unless FORCE_REBUILD).
# ──────────────────────────────────────────────────────────────────────────────
SHAP_PKL_PATH       <- "shap_values_for_sharing.pkl"
SHAP_D_PNG_PATH     <- "SHAP_codonOpt_cutoff_panel.png"
SHAP_D_PDF_PATH     <- "SHAP_codonOpt_cutoff_panel.pdf"
# This standalone SHAP plot is NOT part of the merged figure. It is skipped by
# default so the figure needs no Python. Set BUILD_SHAP_STANDALONE <- TRUE (and
# point SHAP_PYTHON_BIN at a Python with numpy/pandas/matplotlib/statsmodels)
# only if you want the separate SHAP_codonOpt_cutoff_panel.{png,pdf}.
BUILD_SHAP_STANDALONE <- FALSE
FORCE_SHAP_D_REBUILD  <- FALSE

# Optional: explicit path to a Python with numpy/pandas/matplotlib/statsmodels
# installed. NULL = auto-detect via Sys.which("python3"). Same convention as
# the merged Figure 2 script.
SHAP_PYTHON_BIN     <- NULL

# Sizing knobs for the SHAP cutoff figure (wider aspect so it fills its
# full-width row in the merged figure without big margins on the sides).
SHAP_D_FIGSIZE_W    <- 22    # wide so it fills a full row
SHAP_D_FIGSIZE_H    <- 10
SHAP_D_DPI          <- 200
SHAP_D_LABEL_FS     <- 28    # axis labels & title
SHAP_D_TICK_FS      <- 22
SHAP_D_LEGEND_FS    <- 20
SHAP_D_CUT_LABEL_FS <- 22    # green cutoff annotation text

needs_rebuild_D <- BUILD_SHAP_STANDALONE && (FORCE_SHAP_D_REBUILD ||
  !file.exists(SHAP_D_PNG_PATH) ||
  (file.exists(SHAP_PKL_PATH) &&
   file.info(SHAP_PKL_PATH)$mtime > file.info(SHAP_D_PNG_PATH)$mtime))

if (needs_rebuild_D) {
  cat("[Panel D] Regenerating SHAP cutoff scatter from",
      SHAP_PKL_PATH, "...\n")

  if (!file.exists(SHAP_PKL_PATH)) {
    stop(
      "[Panel D] SHAP pickle not found at: ", SHAP_PKL_PATH, "\n",
      "  Working directory is currently: ", getwd(), "\n",
      "  Either move the pickle into the working directory, or edit\n",
      "  SHAP_PKL_PATH at the top of the Panel D block."
    )
  }

  py_script <- tempfile(fileext = ".py")
  writeLines(sprintf('
import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess

PKL_PATH = r"%s"
OUT_PNG  = r"%s"
OUT_PDF  = r"%s"
FIGSIZE  = (%g, %g)
DPI      = %d
LABEL_FS = %g
TICK_FS  = %g
LEGEND_FS = %g
CUT_LABEL_FS = %g

with open(PKL_PATH, "rb") as f:
    bundle = pickle.load(f)
shap_arr = bundle["shap_values"]
X        = bundle["feature_data"]
if isinstance(shap_arr, pd.DataFrame): shap_arr = shap_arr.values

FEATURE = "CodonOptimalityFraction_PTCpm100nt"
fidx = X.columns.get_loc(FEATURE)
x   = X[FEATURE].to_numpy().astype(float)
shp = shap_arr[:, fidx]

def loess_smooth(xv, yv, frac=0.2):
    fit = lowess(yv, xv, frac=frac, return_sorted=True)
    return fit[:, 0], fit[:, 1]
def find_zero_crossings(xs, ys):
    sign = np.sign(ys); sign[sign == 0] = 1
    idx = np.where(np.diff(sign) != 0)[0]
    out = []
    for i in idx:
        x0, x1, y0, y1 = xs[i], xs[i+1], ys[i], ys[i+1]
        if y1 != y0:
            out.append(x0 - y0 * (x1 - x0) / (y1 - y0))
    return np.array(out)
def bootstrap_crossings(xv, yv, frac=0.2, n_boot=300, seed=42):
    rng = np.random.default_rng(seed); n = len(xv); all_c = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        try:
            xs, ys = loess_smooth(xv[idx], yv[idx], frac=frac)
            all_c.extend(find_zero_crossings(xs, ys))
        except Exception:
            continue
    return np.array(all_c)
def cluster(values, gap, min_size=5):
    if len(values) == 0: return []
    s = np.sort(values); clusters = []; current = [s[0]]
    for v in s[1:]:
        if v - current[-1] <= gap: current.append(v)
        else: clusters.append(current); current = [v]
    clusters.append(current)
    return [(np.median(c), np.percentile(c, 2.5), np.percentile(c, 97.5))
            for c in clusters if len(c) >= min_size]

DATA_RANGE = float(x.max() - x.min())
cs = bootstrap_crossings(x, shp, frac=0.2, n_boot=300, seed=42)
clusters_pooled = cluster(cs, gap=DATA_RANGE * 0.05)
print(f"Clusters: {clusters_pooled}")

fig, ax = plt.subplots(figsize=FIGSIZE, facecolor="white", dpi=DPI)
ax.set_facecolor("white")
ax.scatter(x, shp, color="grey", s=24, alpha=0.40,
           edgecolors="white", linewidths=0.3)
xs_full, ys_full = loess_smooth(x, shp, frac=0.2)
ax.plot(xs_full, ys_full, color="#C0392B", linewidth=4.0,
        label="LOESS (frac=0.2)")
ax.axhline(0, color="grey", linestyle="--", linewidth=1.4)

ymin_full, ymax_full = float(shp.min()), float(shp.max())
y_off = [0.96, 0.85]
for i, (med, lo, hi) in enumerate(clusters_pooled):
    ax.axvline(med, color="#27AE60", linewidth=3.0, alpha=0.85)
    y_lab = ymin_full + (ymax_full - ymin_full) * y_off[i %% 2]
    ax.text(med, y_lab,
            f"+{med:.3f}\\n95%% CI [{lo:.2f}, {hi:.2f}]",
            color="#27AE60", fontsize=CUT_LABEL_FS, fontweight="bold",
            ha="center", va="top",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="#27AE60", linewidth=1.4, alpha=0.92))

ax.set_xlabel("Codon optimality (PTC \u00b1 100 nt)",
              fontsize=LABEL_FS, fontweight="bold")
ax.set_ylabel("SHAP value", fontsize=LABEL_FS, fontweight="bold")
ax.set_title("D. Cutoff Analysis: Bootstrap Zero Crossings",
             fontsize=LABEL_FS, fontweight="bold", loc="left")
ax.tick_params(axis="both", labelsize=TICK_FS, width=2.0, length=7,
               colors="#333333")
ax.legend(loc="lower right", fontsize=LEGEND_FS, frameon=True,
          edgecolor="#333333")
ax.grid(True, alpha=0.3, linewidth=0.6)
ax.set_axisbelow(True)
for sp in ["top","right","bottom","left"]:
    ax.spines[sp].set_linewidth(1.8); ax.spines[sp].set_color("#333333")

plt.tight_layout()
fig.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight",
            facecolor="white", edgecolor="none")
fig.savefig(OUT_PDF,            bbox_inches="tight",
            facecolor="white", edgecolor="none")
plt.close(fig)
print(f"Saved: {OUT_PNG}")
', SHAP_PKL_PATH, SHAP_D_PNG_PATH, SHAP_D_PDF_PATH,
   SHAP_D_FIGSIZE_W, SHAP_D_FIGSIZE_H, SHAP_D_DPI,
   SHAP_D_LABEL_FS, SHAP_D_TICK_FS, SHAP_D_LEGEND_FS, SHAP_D_CUT_LABEL_FS),
    py_script)

  if (!is.null(SHAP_PYTHON_BIN) && nzchar(SHAP_PYTHON_BIN)) {
    py_bin <- path.expand(SHAP_PYTHON_BIN)
    if (!file.exists(py_bin))
      stop("[Panel D] SHAP_PYTHON_BIN points to a non-existent file: ", py_bin)
  } else {
    py_bin <- Sys.which("python3")
    if (py_bin == "") py_bin <- Sys.which("python")
    if (py_bin == "")
      stop("[Panel D] No python3/python found on PATH. Set SHAP_PYTHON_BIN.")
  }
  cat("[Panel D] Using Python at:", py_bin, "\n")

  status <- system2(py_bin, args = py_script, stdout = "", stderr = "")
  if (status != 0) {
    stop(
      "[Panel D] SHAP cutoff regeneration failed (exit ", status, ").\n",
      "Most common cause: ", py_bin,
      " is missing numpy/pandas/matplotlib/statsmodels.\n",
      "Fix: ", py_bin, " -m pip install numpy pandas matplotlib statsmodels\n",
      "Or set SHAP_PYTHON_BIN to a Python that has them."
    )
  }
  cat("[Panel D] SHAP cutoff scatter rebuilt.\n")
} else if (BUILD_SHAP_STANDALONE) {
  cat("[Panel D] Reusing existing", SHAP_D_PNG_PATH,
      "(set FORCE_SHAP_D_REBUILD <- TRUE to rebuild).\n")
} else {
  cat("[SHAP standalone] Skipped (not part of this figure).\n")
}

# Embed the PNG as a standalone object. NOTE: this is NOT placed in the merged
# figure (the SHAP codon-optimality panel is omitted); kept only for optional
# standalone use. Renamed from pD to avoid clashing with the re-lettered
# lollipop panels (D/E/F) below.
if (BUILD_SHAP_STANDALONE && file.exists(SHAP_D_PNG_PATH))
p_shap_standalone <- ggdraw() +
  draw_image(SHAP_D_PNG_PATH) +
  draw_label("", x = 0.005, y = 0.99,
             hjust = 0, vjust = 1,
             fontface = "bold",
             size = FONT + 28,
             color = "grey10")

# ══════════════════════════════════════════════════════════════════════════════
# PANELS D / E / F — RBP motif lollipops
#   All RBP statistics computed on the variant-level subset (df_var).
# ══════════════════════════════════════════════════════════════════════════════
make_rbp_lollipop <- function(rbp_col_vec, col_prefix,
                               n_top = 10,
                               bonf_only = FALSE,
                               tag_label, title_label,
                               caption_label = NULL,
                               col_inhibit = "#2980B9",
                               col_promote = "#27AE60") {

  cat(sprintf("Computing RBP statistics for %s (%d columns)...\n",
              col_prefix, length(rbp_col_vec)))

  rbp_stats <- map_dfr(rbp_col_vec, function(col) {
    has <- df_var %>% filter(.data[[col]] == 1) %>% pull(ALLELE.RAT)
    no  <- df_var %>% filter(.data[[col]] == 0) %>% pull(ALLELE.RAT)
    if (length(has) < 10) return(NULL)
    wt <- wilcox.test(has, no, alternative = "two.sided")
    tibble(rbp   = sub("^[^.]+\\.", "", col),
           p_raw = wt$p.value,
           n_has = length(has),
           diff  = mean(has, na.rm = TRUE) - mean(no, na.rm = TRUE))
  }) %>%
    mutate(
      p_bonf     = pmin(p_raw * n(), 1),
      sig_bonf   = case_when(p_bonf < .0001 ~ "****", p_bonf < .001 ~ "***",
                              p_bonf < .01   ~ "**",   p_bonf < .05  ~ "*",
                              TRUE           ~ "ns"),
      sig_raw    = case_when(p_raw  < .0001 ~ "****", p_raw  < .001 ~ "***",
                              p_raw  < .01   ~ "**",   p_raw  < .05  ~ "*",
                              TRUE           ~ "ns")
    )

  if (bonf_only) {
    plot_df <- rbp_stats %>%
      filter(sig_bonf != "ns") %>%
      arrange(diff) %>%
      mutate(direction  = if_else(diff < 0, "Inhibits NMD", "Promotes NMD"),
             show_label = sig_bonf,
             rbp        = factor(rbp, levels = rbp))

    if (nrow(plot_df) == 0) {
      message(sprintf("No Bonferroni-significant motifs for %s", col_prefix))
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5,
                        label = paste0("No Bonferroni-significant RBP motifs\nin ",
                                       title_label),
                        size = FONT / 3, color = "grey40", hjust = 0.5) +
               labs(tag = tag_label, title = title_label) +
               theme_void() +
               theme(plot.tag   = element_text(size = FONT + 20, face = "bold"),
                     plot.title = element_text(size = FONT + 4, face = "bold",
                                               color = TITLE_COL, hjust = 0.5)))
    }

  } else {
    top_neg <- rbp_stats %>% filter(diff < 0) %>%
      arrange(p_bonf) %>% slice_head(n = n_top)
    top_pos <- rbp_stats %>% filter(diff > 0) %>%
      arrange(p_raw)  %>% slice_head(n = n_top)

    plot_df <- bind_rows(
      top_neg %>% mutate(direction = "Inhibits NMD", show_label = sig_bonf),
      top_pos %>% mutate(direction = "Promotes NMD", show_label = sig_raw)
    ) %>%
      arrange(diff) %>%
      mutate(rbp = factor(rbp, levels = rbp))
  }

  x_range   <- max(abs(plot_df$diff)) * 1.45
  x_breaks  <- pretty(c(-x_range, x_range), n = 6)

  ggplot(plot_df, aes(y = rbp)) +
    geom_vline(xintercept = 0, color = "grey25", linewidth = 2.4) +
    annotate("rect", xmin = -Inf, xmax = 0,
             ymin = -Inf, ymax = Inf, fill = "#EBF5FB", alpha = 0.45) +
    annotate("rect", xmin = 0,    xmax = Inf,
             ymin = -Inf, ymax = Inf, fill = "#EAFAF1", alpha = 0.45) +
    geom_segment(aes(x = 0, xend = diff, y = rbp, yend = rbp, color = direction),
                 linewidth = 7.0, alpha = 0.85) +
    geom_point(aes(x = diff, color = direction), size = 22, alpha = 0.95) +
    geom_text(aes(x = diff,
                  label  = paste0("n=", format(n_has, big.mark = ",")),
                  hjust  = if_else(diff < 0, 1.18, -0.18)),
              size = 16, fontface = "italic", color = "grey35") +
    geom_text(aes(x = diff + sign(diff) * 0.0005,
                  label  = show_label,
                  hjust  = if_else(diff < 0, -0.18, 1.18)),
              size = 24, fontface = "bold", color = "grey15") +
    annotate("text",
             x = -x_range * 0.50, y = nrow(plot_df) + 1.4,
             label    = "Inhibits NMD",
             fontface = "bold.italic", size = 22,
             color    = col_inhibit, hjust = 0.5) +
    annotate("text",
             x = x_range * 0.50, y = nrow(plot_df) + 1.4,
             label    = "Promotes NMD",
             fontface = "bold.italic", size = 22,
             color    = col_promote, hjust = 0.5) +
    scale_color_manual(
      values = setNames(c(col_inhibit, col_promote),
                        c("Inhibits NMD", "Promotes NMD")),
      guide  = "none") +
    scale_x_continuous(breaks = x_breaks,
                       limits = c(-x_range, x_range),
                       labels = function(x) sprintf("%.2f", x)) +
    labs(x       = "Mean NMD efficiency difference (motif present \u2212 absent)",
         y       = NULL,
         title   = title_label,
         tag     = tag_label) +
    base_theme +
    coord_cartesian(clip = "off") +
    theme(
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.7),
      panel.grid.major.x = element_line(color = "grey91", linewidth = 0.8),
      axis.text.y        = element_text(size = RBP_MOTIF_SIZE, face = "bold",
                                        color = "grey10"),
      plot.margin        = margin(30, 40, 18, 22)
    )
}

pD <- make_rbp_lollipop(
  rbp_col_vec   = rbp_cols_ejc,
  col_prefix    = "ptc_to_ejc",
  n_top         = 10,
  bonf_only     = FALSE,
  tag_label     = "D",
  title_label   = "Top 10 RBP Motifs in PTC-to-EJC: Inhibiting vs Promoting NMD"
)

pE <- make_rbp_lollipop(
  rbp_col_vec   = rbp_cols_ptc100,
  col_prefix    = "ptcpm100",
  bonf_only     = TRUE,
  tag_label     = "E",
  title_label   = "Bonferroni-Significant RBP Motifs — PTC \u00b1 100 nt"
)

pF <- make_rbp_lollipop(
  rbp_col_vec   = rbp_cols_ejc100,
  col_prefix    = "ejcpm100",
  bonf_only     = TRUE,
  tag_label     = "F",
  title_label   = "Bonferroni-Significant RBP Motifs — EJC \u00b1 100 nt"
)

# ══════════════════════════════════════════════════════════════════════════════
# COMBINE
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
# COMBINE
# The SHAP codon-optimality scatter (built by the block above) is NOT lettered
# in this figure — it is omitted from the merged layout and produced only as a
# standalone file (SHAP_codonOpt_cutoff_panel.{png,pdf}). The six rendered
# panels are therefore lettered consecutively A–F with no gap.
# ══════════════════════════════════════════════════════════════════════════════
combined <- pA / pB / pC / pD / pE / pF +
  plot_layout(heights = c(1.0,   # A — wide hex plot
                          1.0,   # B — 4-facet scatter
                          0.9,   # C — narrow conservation bars
                          1.2,   # D — 20 RBP rows, taller
                          1.0,   # E — 17 RBP rows
                          1.0)) &  # F — 16 RBP rows
  theme(plot.margin = margin(14, 18, 14, 18))

# Canvas height reduced now that Panel D no longer takes a full row.
ggsave("NMD_PTCtoEJC_Features_Figure_revised.pdf", combined, device = pdf_device,
       width = 80, height = 132, dpi = 400, limitsize = FALSE)
# PNG preview at 100 dpi (8,000 x 13,200 px). 300 dpi on this 80 x 132 in
# canvas is ~950 megapixels and exhausts memory; the PDF is the vector master.
ggsave("NMD_PTCtoEJC_Features_Figure_revised.png", combined,
       width = 80, height = 132, dpi = 100, limitsize = FALSE)

print(combined)
cat("Saved: NMD_PTCtoEJC_Features_Figure_revised.pdf / .png\n")
