# ══════════════════════════════════════════════════════════════════════════════
# NMD — Supplemental Figure: 5'UTR Features vs NMD Efficiency  (v4)
#
#   Panel A: 5'UTR AU content tertile         (Low / Medium / High)
#   Panel B: 5'UTR UC content tertile         (Low / Medium / High)
#   Panel C: 5'UTR uORF — All PTC-variants    (Without / With uORF)
#   Panel D: 5'UTR introns                    (Without / With 5'UTR intron)
#   Panel E: 5'UTR PhastCons conservation     (Low / Medium / High)
#
# Changes applied in v4 (vs v3):
#   1) Filter PTC-variants to ALLELE.RAT > 0.35
#   2) GC-content panel removed (kept AU and UC)
#   3) Upstream-PTC uORF panel removed
#   4) No gray italic subtitle/caption labels under titles or panels
#   5) Larger fonts, bolder lines
#   6) Transcript-level analysis = median ALLELE.RAT per transcript
#      (no bootstrapping — one row per transcript)
#
# v5 (corrected dataset, October 2026):
#   - input switched to TOPMed_stopgain_September25_corrected_readyformodel.csv
#   - filter ALLELE.RAT >= 0.35 (was > 0.35), matching all other figures
#   - stops with a clear message if any 5'UTR column is missing, listing the
#     closest column names in the file
#   - prints n per group for every panel (for the legend)
#   - output: SupplementalFigure6_revised.pdf / .png
#
# Required file:
#   TOPMed_stopgain_September25_corrected_readyformodel.csv
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggforce)
library(patchwork)
library(scales)

# Resolve namespace conflicts (defensive)
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
distinct  <- dplyr::distinct
left_join <- dplyr::left_join
arrange   <- dplyr::arrange
count     <- dplyr::count
rename    <- dplyr::rename
slice     <- dplyr::slice
drop_na   <- tidyr::drop_na

set.seed(42)

# ══════════════════════════════════════════════════════════════════════════════
# COLOURS & THEME  (bigger fonts, bolder lines)
# ══════════════════════════════════════════════════════════════════════════════
TITLE_COL <- "#2E6DA4"
FONT      <- 26                         # base font size (was ~17 effectively)

# Tertile palette for AU / UC / PhastCons (3-level)  — picked to mirror image
pal_AU       <- c(Low = "#BBC9D6", Medium = "#7FA4C6", High = "#5C7E9E")
pal_UC       <- c(Low = "#F2D9C2", Medium = "#E5A678", High = "#A1764D")
pal_phast    <- c(Low = "#BFBFBF", Medium = "#7FB4D6", High = "#2F6E94")

# Binary palettes for uORF / intron panels
pal_uorf     <- c(`Without uORF` = "#CFE3F2", `With uORF` = "#C25F5F")
pal_intron   <- c(`Without\n5'UTR intron` = "#CFE3F2",    # v5: two-line labels
                  `With\n5'UTR intron`    = "#A678C2")

base_theme <- theme_classic(base_size = FONT) +
  theme(
    legend.position    = "none",
    axis.text.x        = element_text(size = FONT + 2, face = "bold",
                                      color = "grey10", lineheight = 1.1),
    axis.text.y        = element_text(size = FONT + 1, face = "bold",
                                      color = "grey15"),
    axis.title.x       = element_text(size = FONT + 4, face = "bold",
                                      color = "grey10",
                                      margin = margin(t = 14)),
    axis.title.y       = element_text(size = FONT + 4, face = "bold",
                                      color = "grey10",
                                      margin = margin(r = 14)),
    axis.line          = element_line(color = "grey15", linewidth = 1.20),
    axis.ticks         = element_line(color = "grey20", linewidth = 0.95),
    axis.ticks.length  = unit(5, "pt"),
    panel.grid.major.y = element_line(color = "grey91", linewidth = 0.55),
    panel.grid.major.x = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.tag           = element_text(size = FONT + 16, face = "bold",
                                      color = "grey10"),
    plot.title         = element_text(size = FONT + 6, face = "bold",
                                      color = TITLE_COL, hjust = 0.5,
                                      margin = margin(b = 10)),
    # No subtitle / no caption used in panels (per request)
    plot.subtitle      = element_blank(),
    plot.caption       = element_blank(),
    plot.margin        = margin(18, 22, 18, 22)
  )

ref05 <- geom_hline(yintercept = 0.5, linetype = "dashed",
                    color = "grey35", linewidth = 1.0)

y_main <- scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                             limits = c(-0.04, 1.36), expand = c(0, 0))   # v5: was 1.30

# ══════════════════════════════════════════════════════════════════════════════
# LOAD & PREP — apply ALLELE.RAT > 0.35 filter, then collapse to transcript level
# ══════════════════════════════════════════════════════════════════════════════
cat("Loading data...\n")
df_raw <- read.csv("TOPMed_stopgain_September25_corrected_readyformodel.csv",
                   stringsAsFactors = FALSE)

# (1) Apply ALLELE.RAT >= 0.35 filter to PTC-variants (v5: >=, was >)
df_var <- df_raw %>%
  filter(!is.na(ALLELE.RAT), ALLELE.RAT >= 0.35,
         !is.na(TxName))

cat(sprintf("Variants after ALLELE.RAT >= 0.35 filter: %d (main figures: 5,749)\n",
            nrow(df_var)))

# Helper: collapse variants -> transcripts by median ALLELE.RAT, carrying
# transcript-level features (assumed constant within TxName).
# `feature_cols` = the per-transcript feature column(s) to keep.
collapse_to_tx <- function(d, feature_cols) {
  d %>%
    group_by(TxName, across(all_of(feature_cols))) %>%
    summarise(NMD = median(ALLELE.RAT, na.rm = TRUE),
              n_var = dplyr::n(),
              .groups = "drop")
}

# Tertile helper — Low / Medium / High by 33rd & 67th percentiles
make_tertile <- function(x) {
  q <- quantile(x, probs = c(1/3, 2/3), na.rm = TRUE)
  cut(x,
      breaks         = c(-Inf, q[1], q[2], Inf),
      labels         = c("Low", "Medium", "High"),
      include.lowest = TRUE)
}

# Conservation bin helper (fixed cut-points used in the original figure)
make_phast_bin <- function(x) {
  cut(x,
      breaks         = c(-Inf, 0.1, 0.5, Inf),
      labels         = c("Low", "Medium", "High"),
      include.lowest = TRUE,
      right          = FALSE)
}

# ══════════════════════════════════════════════════════════════════════════════
# Build per-panel transcript-level data frames
# (Adjust the *_col strings below if your CSV uses different column names.)
# ══════════════════════════════════════════════════════════════════════════════

# --- column-name expectations (edit if needed) ---
COL_AU      <- "fiveUTR_AU_content"               # numeric, 5'UTR AU fraction
COL_UC      <- "fiveUTR_UC_content"               # numeric, 5'UTR UC fraction
COL_UORF    <- "fiveutrseqs.uORF"                 # 0/1, T/F, or count → coerced
COL_INTRON  <- "fiveUTR.introns"                  # 0/1, T/F, or count → coerced
COL_PHAST   <- "phastcons_utr5_first200_median"   # mean phastCons, first 200 nt

# Robust coercion of "has feature?" columns to a clean logical.
# Treats NA / 0 / FALSE / "" / "none" as FALSE; everything else as TRUE.
# v5: text values are now read by meaning. Previously any non-empty text
# counted as TRUE, so a value such as "No 5UTR intron" was coded as HAVING an
# intron and every transcript fell into one group (wilcox.test then fails with
# "not enough 'y' observations").
to_has_flag <- function(x) {
  if (is.logical(x))    return(x)
  if (is.numeric(x))    return(x > 0)
  s <- tolower(trimws(as.character(x)))
  neg <- is.na(s) | s %in% c("", "0", "false", "f", "no", "n", "none", "na") |
         grepl("\\bno\\b|\\bnot\\b|without|absent|^none", s)
  num <- suppressWarnings(as.numeric(s))
  ifelse(!is.na(num), num > 0, !neg)
}
# Show exactly how each raw value was coded, so the grouping can be checked
show_flag_map <- function(d, col) {
  cat(sprintf("\n%s -> coded as:\n", col))
  m <- d %>% count(raw = as.character(.data[[col]])) %>%
    mutate(has_feature = to_has_flag(raw))
  print(as.data.frame(m), row.names = FALSE)
}

# Helper: safe presence check
has_col <- function(d, x) x %in% names(d)

# v5: stop early, with suggestions, if a 5'UTR column is missing
req5  <- c(COL_AU, COL_UC, COL_UORF, COL_INTRON, COL_PHAST)
miss5 <- req5[!sapply(req5, has_col, d = df_var)]
if (length(miss5) > 0)
  stop("Missing column(s): ", paste(miss5, collapse = ", "),
       "\nEdit the COL_* settings above. Columns mentioning 5'UTR: ",
       paste(grep("5utr|fiveutr|utr5|uorf", names(df_var), ignore.case = TRUE,
                  value = TRUE), collapse = ", "))

show_flag_map(df_var, COL_UORF)
show_flag_map(df_var, COL_INTRON)

# ---------- Panel A : AU tertile ----------
df_A_var <- df_var %>% filter(!is.na(.data[[COL_AU]]))
df_A_tx  <- collapse_to_tx(df_A_var, COL_AU) %>%
  mutate(grp = make_tertile(.data[[COL_AU]])) %>%
  drop_na(grp)

# ---------- Panel B : UC tertile (was C) ----------
df_B_var <- df_var %>% filter(!is.na(.data[[COL_UC]]))
df_B_tx  <- collapse_to_tx(df_B_var, COL_UC) %>%
  mutate(grp = make_tertile(.data[[COL_UC]])) %>%
  drop_na(grp)

# ---------- Panel C : uORF, all PTC-variants (was D) ----------
df_C_var <- df_var %>%
  # v5: blank = no annotated 5'UTR (92 variants); excluded, not "Without uORF"
  filter(!is.na(.data[[COL_UORF]]), trimws(as.character(.data[[COL_UORF]])) != "") %>%
  mutate(.uorf_flag = to_has_flag(.data[[COL_UORF]]))
df_C_tx  <- collapse_to_tx(df_C_var, ".uorf_flag") %>%
  mutate(grp = factor(ifelse(.uorf_flag, "With uORF", "Without uORF"),
                      levels = c("Without uORF", "With uORF"))) %>%
  drop_na(grp)

# ---------- Panel D : 5'UTR intron (was E, was F in original) ----------
# v5: in the corrected file fiveUTR.introns only records presence
# ("There is a 5UTR intron"); transcripts WITHOUT a 5'UTR intron are blank (NA).
# They must be kept as "Without", not dropped (dropping them left only the
# with-intron group and the Wilcoxon test failed).
df_D_var <- df_var %>%
  mutate(.intron_flag = !is.na(.data[[COL_INTRON]]) & to_has_flag(.data[[COL_INTRON]]))
df_D_tx  <- collapse_to_tx(df_D_var, ".intron_flag") %>%
  mutate(grp = factor(ifelse(.intron_flag,
                             "With\n5'UTR intron", "Without\n5'UTR intron"),
                      levels = c("Without\n5'UTR intron",
                                 "With\n5'UTR intron"))) %>%
  drop_na(grp)

# ---------- Panel E : PhastCons (was F, was G in original) ----------
df_E_var <- df_var %>% filter(!is.na(.data[[COL_PHAST]]))
df_E_tx  <- collapse_to_tx(df_E_var, COL_PHAST) %>%
  mutate(grp = make_phast_bin(.data[[COL_PHAST]])) %>%
  drop_na(grp)

cat("Transcript counts per panel:\n")
cat(sprintf("  A (AU)        : %d transcripts\n", nrow(df_A_tx)))
cat(sprintf("  B (UC)        : %d transcripts\n", nrow(df_B_tx)))
cat(sprintf("  C (uORF all)  : %d transcripts\n", nrow(df_C_tx)))
cat(sprintf("  D (intron)    : %d transcripts\n", nrow(df_D_tx)))
cat(sprintf("  E (phastCons) : %d transcripts\n", nrow(df_E_tx)))

# v5: n per group (for the legend)
cat("\n── n per group (transcripts) ──\n")
for (nm in c("A", "B", "C", "D", "E")) {
  d <- get(paste0("df_", nm, "_tx"))
  cat(nm, ": ", paste(paste0(names(table(d$grp)), " = ", table(d$grp)),
                      collapse = "; "), "\n", sep = "")
}

# ══════════════════════════════════════════════════════════════════════════════
# Plot helpers
# ══════════════════════════════════════════════════════════════════════════════

# Tertile / 3-group violin+box plot (used in A, B)
plot_tertile <- function(d, palette, panel_tag, panel_title, x_lab) {

  comparisons <- list(c("Low","Medium"), c("Medium","High"), c("Low","High"))
  y_brk       <- c(1.05, 1.12, 1.19)   # v5: kept below the n labels

  n_lab <- d %>% count(grp) %>%
    mutate(label = format(n, big.mark = ","))

  st <- d %>%
    wilcox_test(NMD ~ grp, comparisons = comparisons) %>%
    add_significance("p") %>%
    add_xy_position(x = "grp") %>%
    mutate(y.position = y_brk, grp = group1)   # v5: inherited aes

  ggplot(d, aes(x = grp, y = NMD, fill = grp, color = grp)) +
    geom_violin(alpha = 0.30, width = 0.85, trim = TRUE, linewidth = 0.6) +
    geom_boxplot(width = 0.32, outlier.shape = NA, linewidth = 0.95,
                 fill = "white", alpha = 0.95, color = "grey10") +
    stat_summary(fun = median, geom = "point",
                 size = 5.2, shape = 21,
                 fill = "white", color = "black", stroke = 1.4) +
    ref05 +
    geom_text(data = n_lab,
              aes(x = grp, y = 1.31, label = label),   # v5: was 1.27
              inherit.aes = FALSE, size = 6.0,
              fontface = "bold", color = "grey10") +
    stat_pvalue_manual(st, label = "p.signif",
                       tip.length = 0.012, bracket.size = 0.85,
                       size = 7.0, color = "grey15") +
    scale_fill_manual(values  = palette) +
    scale_color_manual(values = palette) +
    y_main +
    labs(x = x_lab, y = "NMD efficiency",
         tag = panel_tag, title = panel_title) +
    base_theme
}

# Binary violin+box plot (used in C, D, E)
plot_binary <- function(d, palette, panel_tag, panel_title, x_lab) {

  lvls        <- levels(d$grp)
  comparisons <- list(lvls)
  n_lab <- d %>% count(grp) %>%
    mutate(label = format(n, big.mark = ","))

  # v5: skip the test (with a message) if either group has < 2 transcripts
  if (nrow(n_lab) < 2 || any(n_lab$n < 2)) {
    message(sprintf("[Panel %s] a group has < 2 transcripts (%s); test skipped.",
                    panel_tag, paste(n_lab$grp, n_lab$n, sep = " = ", collapse = ", ")))
    st <- NULL
  } else {
    st <- d %>%
      wilcox_test(NMD ~ grp, comparisons = comparisons) %>%
      add_significance("p") %>%
      add_xy_position(x = "grp") %>%
      mutate(y.position = 1.10, grp = group1)   # v5: inherited aes
  }

  ggplot(d, aes(x = grp, y = NMD, fill = grp, color = grp)) +
    geom_violin(alpha = 0.30, width = 0.85, trim = TRUE, linewidth = 0.6) +
    geom_boxplot(width = 0.30, outlier.shape = NA, linewidth = 0.95,
                 fill = "white", alpha = 0.95, color = "grey10") +
    stat_summary(fun = median, geom = "point",
                 size = 5.2, shape = 21,
                 fill = "white", color = "black", stroke = 1.4) +
    ref05 +
    geom_text(data = n_lab,
              aes(x = grp, y = 1.21, label = label),
              inherit.aes = FALSE, size = 6.0,
              fontface = "bold", color = "grey10") +
    { if (!is.null(st)) stat_pvalue_manual(st, label = "p.signif",
                       tip.length = 0.012, bracket.size = 0.85,
                       size = 7.0, color = "grey15") } +
    scale_fill_manual(values  = palette) +
    scale_color_manual(values = palette) +
    y_main +
    labs(x = x_lab, y = "NMD efficiency",
         tag = panel_tag, title = panel_title) +
    base_theme
}

# Bar plot for PhastCons (Panel F) — mean ± 95% CI
plot_phast_bar <- function(d, palette, panel_tag, panel_title, x_lab) {

  comparisons <- list(c("Low","Medium"), c("Medium","High"), c("Low","High"))
  y_brk       <- c(0.86, 0.92, 0.98)

  summ <- d %>%
    group_by(grp) %>%
    summarise(mean_nmd = mean(NMD, na.rm = TRUE),
              sd_nmd   = sd(NMD,   na.rm = TRUE),
              n        = dplyr::n(),
              se       = sd_nmd / sqrt(n),
              ci       = qt(0.975, df = pmax(n - 1, 1)) * se,
              .groups  = "drop")

  st <- d %>%
    wilcox_test(NMD ~ grp, comparisons = comparisons) %>%
    add_significance("p") %>%
    add_xy_position(x = "grp") %>%
    mutate(y.position = y_brk, grp = group1)   # v5: inherited aes

  ggplot(summ, aes(x = grp, y = mean_nmd, fill = grp)) +
    geom_col(width = 0.72, color = "grey10", linewidth = 0.9) +
    geom_errorbar(aes(ymin = mean_nmd - ci, ymax = mean_nmd + ci),
                  width = 0.22, linewidth = 1.0, color = "grey10") +
    geom_text(aes(y = mean_nmd + ci + 0.04,
                  label = sprintf("%.3f", mean_nmd)),
              size = 6.0, fontface = "bold", color = "grey10") +
    geom_text(aes(y = 0.04, label = format(n, big.mark = ",")),
              size = 5.5, fontface = "bold", color = "white") +
    ref05 +
    stat_pvalue_manual(st, label = "p.signif",
                       tip.length = 0.012, bracket.size = 0.85,
                       size = 7.0, color = "grey15") +
    scale_fill_manual(values = palette) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                       limits = c(0, 1.10), expand = c(0, 0)) +
    labs(x = x_lab, y = "Mean NMD efficiency",
         tag = panel_tag, title = panel_title) +
    base_theme
}

# ══════════════════════════════════════════════════════════════════════════════
# Build panels
# ══════════════════════════════════════════════════════════════════════════════
pA <- plot_tertile(df_A_tx, pal_AU,    "A",
                   "5'UTR AU Content", "5'UTR AU content tertile")

pB <- plot_tertile(df_B_tx, pal_UC,    "B",
                   "5'UTR UC Content", "5'UTR UC content tertile")

pC <- plot_binary (df_C_tx, pal_uorf,  "C",
                   "5'UTR uORF \u2014 All PTC-variants", "5'UTR uORF status")

pD <- plot_binary (df_D_tx, pal_intron,"D",
                   "5'UTR Introns", "5'UTR intron status")

pE <- plot_phast_bar(df_E_tx, pal_phast, "E",
                     "5'UTR PhastCons Conservation",
                     "PhastCons conservation (5'UTR first 200 nt)")

# ══════════════════════════════════════════════════════════════════════════════
# Compose figure: 5 panels
#   Row 1: A | B | C       (3 panels)
#   Row 2:   D | E         (2 panels, centered with side spacers)
# ══════════════════════════════════════════════════════════════════════════════
top_row    <- pA | pB | pC
bottom_row <- plot_spacer() | pD | pE | plot_spacer()

combined <- top_row / bottom_row +
  plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(14, 18, 14, 18))

combined <- combined +
  plot_annotation(
    title = "5'UTR Features vs NMD Efficiency",
    theme = theme(
      plot.title = element_text(size = FONT + 10, face = "bold",
                                color = TITLE_COL, hjust = 0.5,
                                margin = margin(b = 6, t = 4))
    )
  )

ggsave("SupplementalFigure6_revised.pdf", combined,
       width = 30, height = 18, dpi = 300, limitsize = FALSE)
ggsave("SupplementalFigure6_revised.png", combined,
       width = 30, height = 18, dpi = 300, limitsize = FALSE)

cat("Saved: SupplementalFigure6_revised.pdf / .png\n")
