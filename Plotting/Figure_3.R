# ══════════════════════════════════════════════════════════════════════════════
# NMD — Merged Figure: Figure 3 A+B + Penultimate Exon PTC-to-EJC Analysis
# v8: Panel E keeps all penultimate variants (n = 585; subtitle states n);
#     open-ended top quintile labelled "129+"; quintile labels horizontal; value
#     labels / brackets / y-limit placed from the data so nothing overlaps;
#     Panel D marker relabelled "Local minimum" (was "Inflection").
# v7: Panels A, B, C unified into one plot style (violin + box + median)
# (FILTERED VERSION v5 — restricted to variants with ALLELE.RAT >= 0.35,
#                       penultimate set further restricted to transcripts with
#                       >= 3 exons, exon-length groups COLLAPSED to 2 tiers
#                       (0-200 bp vs >200 bp), per-group quintile bars (G row)
#                       REMOVED, global font size and line weights increased,
#                       inflection point on Panel D made more prominent)
#
# Layout:
#   Row 1: A (Exon position) | B (EJC downstream) | C (EJC overlap)
#          — all three: violin + box + median, identical styling
#   Row 2: D (Scatter + loess + inflection point, all penultimate, x<=250)
#           | E (Quintile bar, all penultimate)
#   Row 3: F (0-200 bp scatter) | G (>200 bp scatter)
#   D, F, G share one scatter style
#
# Required file:
#   TOPMed_stopgain_September25_corrected_readyformodel.csv
#
# Filters applied:
#   1. ALLELE.RAT >= 0.35 (excludes low-expression / strongly NMD-triggered
#      measurements that may be noise-dominated)
#   2. Penultimate subset (df_penu): only transcripts with exon_count >= 3
#      (penultimate is undefined / meaningless for transcripts with <3 exons)
#
# Penultimate-exon-length groups (length.mutated.exon, in bp):
#   0-200: penultimate exon up to 200 bp
#   >200:  penultimate exon longer than 200 bp
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
# library(ggdist)    # v7: no longer needed (Panel A is now violin + box)
# library(gghalves)   # no longer needed (see left-side point layer)
library(patchwork)

set.seed(42)

# ══════════════════════════════════════════════════════════════════════════════
# HELPER — local minima from loess curve
# ══════════════════════════════════════════════════════════════════════════════
find_local_minima <- function(x, y, span = 0.75, n_seq = 400) {
  lo    <- loess(y ~ x, span = span)
  x_seq <- seq(min(x), max(x), length.out = n_seq)
  y_hat <- predict(lo, newdata = data.frame(x = x_seq))
  is_min <- c(FALSE,
              y_hat[-c(1, n_seq)] < y_hat[-c(n_seq - 1, n_seq)] &
              y_hat[-c(1, n_seq)] < y_hat[-c(1, 2)],
              FALSE)
  tibble(x_min = x_seq[is_min], y_min = y_hat[is_min])
}

# ══════════════════════════════════════════════════════════════════════════════
# COLOURS & THEME
# ══════════════════════════════════════════════════════════════════════════════
TITLE_COL  <- "#2E6DA4"
FONT       <- 48              # was 36 — bumped for bolder appearance
MIN_COL    <- "#C0392B"
LW_MULT    <- 1.6             # global linewidth multiplier (used below)

# Fig3 A colours
col_blue   <- "#6EB4E8"
col_grey   <- "#909090"
col_teal   <- "#3EC9A7"
col_orange <- "#E8924A"
pal_ejc    <- c("No" = col_orange, "Yes" = col_teal)

# Penultimate exon group colours
col_short     <- "#E67E22"
col_medium    <- "#2980B9"
col_long      <- "#27AE60"
col_verylong  <- "#8E44AD"   # purple — kept for back-compat / unused in v4
pal_exon      <- c("0\u2013200" = col_medium,
                   ">200"       = col_long)

base_theme <- theme_classic(base_size = FONT) +
  theme(
    legend.position    = "none",
    axis.text.x        = element_text(size = FONT + 2, face = "bold",
                                      color = "grey10", lineheight = 1.1),
    axis.text.y        = element_text(size = FONT + 1, face = "bold",
                                      color = "grey15"),
    axis.title.y       = element_text(size = FONT + 3, face = "bold",
                                      margin = margin(r = 14)),
    axis.title.x       = element_blank(),
    axis.line          = element_line(color = "grey20", linewidth = 1.36),
    axis.ticks         = element_line(color = "grey35", linewidth = 0.88),
    panel.grid.major.y = element_line(color = "grey91", linewidth = 0.45),
    panel.grid.major.x = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    strip.text         = element_text(size = FONT + 1, face = "bold",
                                      color = "grey10",
                                      margin = margin(8, 6, 8, 6)),
    strip.background   = element_rect(fill = "grey94", color = "grey60",
                                      linewidth = 1.04),
    plot.tag           = element_text(size = FONT + 28, face = "bold",
                                      color = "grey10"),
    plot.title         = element_text(size = FONT + 3, face = "bold",
                                      color = TITLE_COL, hjust = 0.5,
                                      margin = margin(b = 10)),
    plot.subtitle      = element_text(size = FONT,     face = "bold.italic",
                                      color = "grey40", hjust = 0.5,
                                      margin = margin(b = 8)),
    plot.caption       = element_text(size = FONT - 4,  face = "italic",
                                      color = "grey50",
                                      hjust = 0.5, margin = margin(t = 8)),
    plot.margin        = margin(18, 22, 18, 22),
    panel.spacing      = unit(1.0, "lines")
  )

scatter_theme <- base_theme +
  theme(
    axis.title.x     = element_text(size = FONT + 3, face = "bold",
                                    margin = margin(t = 12)),
    axis.title.y     = element_text(size = FONT + 3, face = "bold",
                                    margin = margin(r = 14)),
    panel.grid.major = element_line(color = "grey91", linewidth = 0.45),
    panel.grid.major.x = element_line(color = "grey91", linewidth = 0.45)
  )

y_nmd_dist <- scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                                   limits = c(-0.04, 1.52), expand = c(0, 0))
y_nmd_scat <- scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                                   limits = c(0, 1.05), expand = c(0, 0))
ref05 <- geom_hline(yintercept = 0.5, linetype = "dashed",
                     color = "grey40", linewidth = 0.96)

# ══════════════════════════════════════════════════════════════════════════════
# LOAD & PREP
# ══════════════════════════════════════════════════════════════════════════════
df_raw <- read.csv("TOPMed_stopgain_September25_corrected_readyformodel.csv",
                   stringsAsFactors = FALSE) %>%
  filter(!is.na(ALLELE.RAT),
         ALLELE.RAT >= 0.35) %>%   # <<< NEW: restrict to variants with allele.rat >= 0.35
  mutate(
    # Panel A — exon position
    exon_pos = case_when(
      last.exon == "lastexon" ~ "Last exon",
      last.exon != "lastexon" &
        penultimate.last50bp == "penultimate.last50bp" ~ "Penultimate",
      TRUE ~ "Upstream"
    ),
    exon_pos = factor(exon_pos,
                       levels = c("Last exon","Penultimate","Upstream")),

    # Panel B — EJC
    ejc = case_when(
      downstream == "downstream of last EJC" ~ "No",
      downstream == "upstream of last EJC"   ~ "Yes"
    ),
    ejc = factor(ejc, levels = c("No","Yes")),

    # Penultimate exon length groups — 2-tier scheme:
    #   0-200      penultimate exon up to 200 bp
    #   >200       penultimate exon longer than 200 bp
    exon_len_grp = case_when(
      penultimate.full == "penultimate" &
        length.mutated.exon <= 200 ~ "0\u2013200",
      penultimate.full == "penultimate" ~ ">200",
      TRUE ~ NA_character_
    ),
    exon_len_grp = factor(exon_len_grp,
                           levels = c("0\u2013200", ">200"))
  )

# Penultimate exon subset — additionally require exon_count >= 3
# (penultimate is undefined for transcripts with fewer than 3 exons)
df_penu <- df_raw %>%
  filter(penultimate.full == "penultimate",
         exon_count >= 3,
         !is.na(PTC.2.EJC), !is.na(length.mutated.exon))

cat(sprintf("All variants: %d\n", nrow(df_raw)))
cat(sprintf("Penultimate:  %d  |  Unique tx: %d\n",
            nrow(df_penu), n_distinct(df_penu$TxName)))

# Filter-impact diagnostics
cat(sprintf("\n[Filter applied: ALLELE.RAT >= 0.35]\n"))
cat("By exon position (filtered set):\n")
print(df_raw %>% count(exon_pos))
cat("By EJC group (filtered set):\n")
print(df_raw %>% count(ejc))
cat("By exon-length group (penultimate, filtered):\n")
print(df_penu %>% count(exon_len_grp))

# ══════════════════════════════════════════════════════════════════════════════
# PANELS A–C — ONE SHARED STYLE (v7, Reviewer request)
#   All three categorical comparisons now use the identical plot type:
#   violin (width-scaled) + narrow boxplot + white median point,
#   same y-axis, same n-labels ("n = …"), same Wilcoxon bracket style,
#   no separately drawn outliers (the violin already shows the full spread).
#   Panel widths are proportional to the number of groups (3 : 2 : 2), so every
#   violin has the same physical width across A, B and C.
# ══════════════════════════════════════════════════════════════════════════════

# Column used for Panel C
df_raw <- df_raw %>%
  mutate(
    ejc_overlap = case_when(
      has_ejc_overlap == "yes" ~ "Yes",
      has_ejc_overlap == "no"  ~ "No",
      TRUE ~ NA_character_
    ),
    ejc_overlap = factor(ejc_overlap, levels = c("No", "Yes"))
  )

pal_A      <- c("Last exon"   = col_blue,
                "Penultimate" = col_grey,
                "Upstream"    = col_teal)
pal_ejc_ov <- c("No" = "#2ECC71", "Yes" = "#9B59B6")

make_violin_panel <- function(data, xvar, pal, comparisons, bracket_y,
                              tag, title) {
  d <- data %>% drop_na(all_of(xvar))
  f <- as.formula(paste("ALLELE.RAT ~", xvar))

  n_lab <- d %>% count(.data[[xvar]]) %>%
    mutate(label = paste0("n = ", format(n, big.mark = ",")))

  st <- d %>%
    wilcox_test(f, comparisons = comparisons) %>%
    add_significance("p") %>%
    add_xy_position(x = xvar) %>%
    mutate(y.position = bracket_y,
           !!xvar := group1)

  ggplot(d, aes(x = .data[[xvar]], y = ALLELE.RAT)) +
    geom_violin(aes(fill = .data[[xvar]], color = .data[[xvar]]),
                alpha = 0.42, width = 0.80, scale = "width",
                trim = TRUE, linewidth = 1.2) +
    geom_boxplot(width = 0.20, outlier.shape = NA, linewidth = 1.6,
                 fill = "white", alpha = 0.88, color = "grey10") +
    stat_summary(fun = median, geom = "point", size = 4.8,
                 shape = 21, fill = "white", color = "black", stroke = 1.3) +
    ref05 +
    geom_text(data = n_lab, aes(x = .data[[xvar]], y = 1.42, label = label),
              inherit.aes = FALSE, size = 14,
              fontface = "bold", color = "grey15") +
    stat_pvalue_manual(st, label = "p.signif",
                       tip.length = 0.012, bracket.size = 1.4,
                       size = 18, color = "grey20") +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) +
    y_nmd_dist +
    labs(y = "NMD efficiency", tag = tag, title = title) +
    base_theme +
    # Same x-label style in A, B and C: slightly smaller and tilted so the
    # category names don't overlap.
    theme(axis.text.x = element_text(size = FONT - 6, face = "bold",
                                     color = "grey10",
                                     angle = 25, hjust = 1, vjust = 1))
}

# PANEL A — PTC-bearing exon position
pA <- make_violin_panel(
  df_raw, "exon_pos", pal_A,
  comparisons = list(c("Last exon", "Penultimate"),
                     c("Penultimate", "Upstream"),
                     c("Last exon", "Upstream")),
  bracket_y = c(1.08, 1.18, 1.28),
  tag = "A", title = "Exon position")

# PANEL B — EJC downstream of PTC (No / Yes)
pB <- make_violin_panel(
  df_raw, "ejc", pal_ejc,
  comparisons = list(c("No", "Yes")), bracket_y = 1.18,
  tag = "B", title = "EJC downstream of PTC")

# PANEL C — EJC overlap (No / Yes)
pC_ejc <- make_violin_panel(
  df_raw, "ejc_overlap", pal_ejc_ov,
  comparisons = list(c("No", "Yes")), bracket_y = 1.18,
  tag = "C", title = "EJC overlap")

# ══════════════════════════════════════════════════════════════════════════════
# PANELS D, F, G — ONE SHARED SCATTER STYLE (v7)
#   Style chosen with SCATTER_STYLE below; every panel gets the same black loess
#   fit (grey 95% CI), dashed 0.5 line, "n = …" top-left, Spearman rho
#   bottom-right. D additionally marks the loess inflection point.
#     D: all penultimate exon variants,    PTC-to-EJC <= 250 nt
#     F: penultimate exon 0-200 bp,        PTC-to-EJC <= 200 nt
#     G: penultimate exon  >200 bp,        PTC-to-EJC <= 400 nt
# ══════════════════════════════════════════════════════════════════════════════
# ── Scatter style for D, F, G ────────────────────────────────────────────────
#   Pick ONE style; D, F and G always share it.
#     "points"  : filled circles + loess (95% CI)
#     "hexbin"  : hexagonal density of variants + loess
#     "density" : 2-D density contours + faint points + loess
#     "binned"  : faint points + binned mean ± 95% CI + loess
#     "pointdensity" : points coloured by local density + loess
#     "quantile": binned median with 25–75% / 10–90% bands + loess
#     "boxbin"  : violin + box per distance bin (matches A–C) + loess
SCATTER_STYLE <- "boxbin"
N_BINS        <- 10      # bins for "binned" style

col_all <- "#34495E"   # slate — all penultimate variants (Panel D)

make_scatter_panel <- function(sub, pt_col, xmax, xbrk, title, tag,
                               span = 0.75, show_inflection = FALSE,
                               style = SCATTER_STYLE) {
  if (nrow(sub) < 5) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = sprintf("%s\n(n = %d, too few)", title, nrow(sub)),
                      size = 8, color = "grey30") +
             labs(tag = tag) + theme_void())
  }
  rho <- cor(sub$PTC.2.EJC, sub$ALLELE.RAT, use = "complete.obs",
             method = "spearman")
  p <- ggplot(sub, aes(x = PTC.2.EJC, y = ALLELE.RAT))

  # ---- data layer (the only part that differs between styles) ----
  if (style == "points") {
    p <- p + geom_point(shape = 21, size = 6.5, fill = pt_col, color = "white",
                        stroke = 0.6, alpha = 0.70)
  } else if (style == "hexbin") {
    p <- p + geom_hex(binwidth = c(xmax / 20, 0.05), color = "white", linewidth = 0.4) +
      scale_fill_gradient(low = alpha(pt_col, 0.15), high = pt_col,
                          name = "Variants",
                          guide = guide_colorbar(barwidth = 1.8, barheight = 22))
  } else if (style == "density") {
    p <- p +
      geom_point(color = pt_col, size = 3, alpha = 0.25, shape = 16) +
      stat_density_2d(aes(fill = after_stat(level)), geom = "polygon",
                      alpha = 0.45, color = pt_col, linewidth = 0.6,
                      contour_var = "ndensity") +
      scale_fill_gradient(low = alpha(pt_col, 0.05), high = pt_col, guide = "none")
  } else if (style == "binned") {
    brks <- seq(0, xmax, length.out = N_BINS + 1)
    bins <- sub %>%
      mutate(bin = cut(PTC.2.EJC, brks, include.lowest = TRUE)) %>%
      drop_na(bin) %>% group_by(bin) %>%
      summarise(x = mean(PTC.2.EJC), m = mean(ALLELE.RAT),
                se = sd(ALLELE.RAT) / sqrt(n()), n = n(), .groups = "drop") %>%
      filter(n >= 3)
    p <- p +
      geom_point(color = pt_col, size = 3.5, alpha = 0.22, shape = 16) +
      geom_errorbar(data = bins, aes(x = x, ymin = m - 1.96 * se, ymax = m + 1.96 * se),
                    inherit.aes = FALSE, width = xmax * 0.018,
                    linewidth = 1.6, color = pt_col) +
      geom_point(data = bins, aes(x = x, y = m), inherit.aes = FALSE,
                 shape = 21, size = 8, fill = pt_col, color = "white", stroke = 1.4)
  } else if (style == "pointdensity") {
    # points coloured by local 2-D density (dense regions darker)
    kd  <- MASS::kde2d(sub$PTC.2.EJC, sub$ALLELE.RAT, n = 120)
    sub$dens <- kd$z[cbind(findInterval(sub$PTC.2.EJC, kd$x),
                           findInterval(sub$ALLELE.RAT, kd$y))]
    sub <- sub %>% arrange(dens)
    p <- ggplot(sub, aes(x = PTC.2.EJC, y = ALLELE.RAT)) +
      geom_point(aes(color = dens), size = 6, alpha = 0.9, shape = 16) +
      scale_color_gradientn(colours = c("grey88", alpha(pt_col, 0.55), pt_col),
                            guide = "none")
  } else if (style %in% c("quantile", "boxbin")) {
    brks <- seq(0, xmax, length.out = N_BINS + 1)
    bw   <- brks[2] - brks[1]
    sub_b <- sub %>%
      mutate(bin = cut(PTC.2.EJC, brks, include.lowest = TRUE),
             bin_mid = brks[as.integer(bin)] + bw / 2) %>%
      drop_na(bin) %>% group_by(bin) %>% filter(n() >= 5) %>% ungroup()
    if (style == "quantile") {
      # median with 25–75% and 10–90% bands per distance bin
      q <- sub_b %>% group_by(bin_mid) %>%
        summarise(q10 = quantile(ALLELE.RAT, .10), q25 = quantile(ALLELE.RAT, .25),
                  med = median(ALLELE.RAT),
                  q75 = quantile(ALLELE.RAT, .75), q90 = quantile(ALLELE.RAT, .90),
                  .groups = "drop")
      p <- p +
        geom_ribbon(data = q, aes(x = bin_mid, ymin = q10, ymax = q90),
                    inherit.aes = FALSE, fill = pt_col, alpha = 0.15) +
        geom_ribbon(data = q, aes(x = bin_mid, ymin = q25, ymax = q75),
                    inherit.aes = FALSE, fill = pt_col, alpha = 0.35) +
        geom_line(data = q, aes(x = bin_mid, y = med), inherit.aes = FALSE,
                  color = pt_col, linewidth = 2.2) +
        geom_point(data = q, aes(x = bin_mid, y = med), inherit.aes = FALSE,
                   shape = 21, size = 7, fill = "white", color = pt_col, stroke = 2)
    } else {
      # one violin + box per distance bin — same look as Panels A–C
      p <- p +
        geom_violin(data = sub_b, aes(x = bin_mid, y = ALLELE.RAT, group = bin_mid),
                    inherit.aes = FALSE, width = bw * 0.92, scale = "width",
                    trim = TRUE, fill = pt_col, alpha = 0.42,
                    color = pt_col, linewidth = 1.0) +
        geom_boxplot(data = sub_b, aes(x = bin_mid, y = ALLELE.RAT, group = bin_mid),
                     inherit.aes = FALSE, width = bw * 0.30, outlier.shape = NA,
                     fill = "white", alpha = 0.88, color = "grey10", linewidth = 1.3) +
        stat_summary(data = sub_b, aes(x = bin_mid, y = ALLELE.RAT, group = bin_mid),
                     inherit.aes = FALSE, fun = median, geom = "point",
                     shape = 21, size = 4.2, fill = "white", color = "black", stroke = 1.2)
    }
  } else stop("Unknown style: ", style)

  # ---- shared layers ----
  p <- p + geom_smooth(method = "loess", formula = y ~ x, span = span,
                       color = "black", linewidth = 2.6, fill = "grey60",
                       alpha = 0.22, se = TRUE)

  if (show_inflection) {
    mins <- find_local_minima(sub$PTC.2.EJC, sub$ALLELE.RAT, span = span)
    if (nrow(mins) > 0) p <- p +
      geom_segment(data = mins, aes(x = x_min, xend = x_min, y = 0, yend = y_min),
                   color = MIN_COL, linetype = "dashed", linewidth = 1.7,
                   inherit.aes = FALSE) +
      geom_point(data = mins, aes(x = x_min, y = y_min),
                 shape = 25, size = 11, fill = MIN_COL, color = "white",
                 stroke = 1.6, inherit.aes = FALSE) +
      geom_label(data = mins,
                 aes(x = x_min + xmax * 0.016, y = y_min - 0.13,
                     # v8: "Local minimum" (it is the loess minimum, not an
                     # inflection point; matches the revised text/response)
                     label = sprintf("Local minimum\n(%.0f nt)", x_min)),
                 hjust = 0, size = 16, fontface = "bold.italic",
                 color = MIN_COL, fill = alpha("white", 0.92),
                 label.size = 0.8, inherit.aes = FALSE)
  }

  p +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "grey40", linewidth = 1.04) +
    annotate("text", x = xmax * 0.02, y = 0.98,
             label = paste0("n = ", format(nrow(sub), big.mark = ",")),
             hjust = 0, size = 18, fontface = "bold", color = "grey20") +
    annotate("text", x = xmax * 0.98, y = 0.08,
             label = sprintf("bolditalic(rho == '%.3f')", rho), parse = TRUE,   # Spearman rho
             hjust = 1, size = 18, fontface = "bold.italic", color = "grey20") +
    # coord_cartesian zooms without dropping data, so hexagons / contours at
    # the plot edges are drawn instead of removed
    scale_x_continuous(breaks = xbrk, expand = c(0.01, 0)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), expand = c(0, 0)) +
    coord_cartesian(xlim = c(-1, xmax * 1.03), ylim = c(0, 1.05)) +
    labs(x = "PTC-to-EJC distance (nt)", y = "NMD efficiency",
         title = title, tag = tag) +
    scatter_theme +
    theme(legend.position = if (style == "hexbin") "right" else "none",
          legend.title = element_text(size = FONT - 4, face = "bold"),
          legend.text  = element_text(size = FONT - 6))
}

# PANEL D — all penultimate exon variants (x <= 250 nt), with inflection point
X_MAX_D <- 250
df_D    <- df_penu %>% filter(PTC.2.EJC <= X_MAX_D)
cat(sprintf("\nPanel D: n=%d (removed %d >%d nt)\n",
            nrow(df_D), sum(df_penu$PTC.2.EJC > X_MAX_D), X_MAX_D))

pD <- make_scatter_panel(df_D, col_all, X_MAX_D, seq(0, X_MAX_D, 50),
                         title = "All penultimate exons", tag = "D",
                         span = 0.75, show_inflection = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# PANEL E — Quintile bar
#   Uses ALL penultimate exon variants (df_penu, no distance cap; n = 585),
#   unlike Panel D, which is windowed to PTC-to-EJC <= X_MAX_D (n = 545).
#   The subtitle states the n so the difference from D is explicit.
#   Set E_MATCH_D <- TRUE to restrict E to Panel D's variants instead.
# ══════════════════════════════════════════════════════════════════════════════
E_MATCH_D <- FALSE   # FALSE = all penultimate variants (585)
df_E <- if (E_MATCH_D) df_D else df_penu

# Quintile breaks computed once and shared by the bars and the Wilcoxon tests
q_brks_E <- quantile(df_E$PTC.2.EJC, probs = seq(0, 1, 0.2), na.rm = TRUE)
if (anyDuplicated(q_brks_E))
  stop("Panel E: tied quintile breaks (", paste(q_brks_E, collapse = ", "),
       "). Too many identical PTC-to-EJC values to form 5 distinct quintiles.")
df_E <- df_E %>%
  mutate(dist_q = cut(PTC.2.EJC, breaks = q_brks_E,
                      include.lowest = TRUE, labels = FALSE)) %>%
  drop_na(dist_q)

stopifnot(!E_MATCH_D || nrow(df_E) == nrow(df_D))
cat(sprintf("Panel E: n=%d (%s)\n", nrow(df_E),
            if (E_MATCH_D) sprintf("same variants as Panel D, <=%d nt", X_MAX_D)
            else "all penultimate, no distance cap"))

bin_E <- df_E %>%
  group_by(dist_q) %>%
  summarise(
    mean_NMD = mean(ALLELE.RAT, na.rm = TRUE),
    se       = sd(ALLELE.RAT,   na.rm = TRUE) / sqrt(n()),
    ci_lo    = mean_NMD - 1.96 * se,
    ci_hi    = mean_NMD + 1.96 * se,
    n        = n(),
    d_min    = min(PTC.2.EJC), d_max = max(PTC.2.EJC),
    d_med    = median(PTC.2.EJC),
    .groups  = "drop"
  ) %>%
  # Range without "nt" (unit is in the axis title). The top quintile is
  # open-ended when uncapped (e.g. 129-4291), which is too wide for a bar
  # label, so it is written as "129+" (its maximum is printed to the console
  # for the figure legend).
  mutate(rng = if_else(dist_q == max(dist_q) & !E_MATCH_D,
                       sprintf("%d+", round(d_min)),
                       sprintf("%d\u2013%d", round(d_min), round(d_max))),
         lab = sprintf("Q%d\n%s\nn = %d", dist_q, rng, n))

cat("Panel E quintiles:\n"); print(bin_E %>% select(dist_q, n, d_min, d_max, mean_NMD))

# Vertical layout derived from the data, so value labels, brackets and the
# panel top never collide:
#   value label  : just above each CI whisker
#   brackets     : two alternating tiers above the tallest value label
TOP_E      <- max(bin_E$ci_hi)
VAL_GAP_E  <- 0.02
BRK_LOW_E  <- TOP_E + 0.10
BRK_HIGH_E <- TOP_E + 0.18
Y_MAX_E    <- TOP_E + 0.29

# Sequential Wilcoxon tests between adjacent quintiles (same quintiles as bars)
stat_E_bar <- df_E %>%
  mutate(dist_q_fac = factor(dist_q)) %>%
  wilcox_test(ALLELE.RAT ~ dist_q_fac,
              comparisons = list(c("1","2"),c("2","3"),c("3","4"),c("4","5"))) %>%
  add_significance("p") %>%
  add_xy_position(x = "dist_q_fac") %>%
  mutate(y.position = c(BRK_LOW_E, BRK_HIGH_E, BRK_LOW_E, BRK_HIGH_E))

# Map stat positions back to bar x-axis
stat_E_bar$x <- stat_E_bar$xmin
stat_E_bar$xmin_bar <- stat_E_bar$xmin
stat_E_bar$xmax_bar <- stat_E_bar$xmax

pE <- ggplot(bin_E,
             aes(x = reorder(lab, dist_q), y = mean_NMD, fill = d_med)) +

  geom_col(alpha = 0.80, color = "grey20", linewidth = 0.8, width = 0.65) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.22, linewidth = 1.76, color = "grey15") +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "grey40", linewidth = 1.04) +
  geom_text(aes(label = sprintf("%.3f", mean_NMD), y = ci_hi + VAL_GAP_E),
            vjust = 0, size = 13, fontface = "bold", color = "grey20") +

  # Sequential p-values between adjacent quintiles
  geom_bracket(data = stat_E_bar,
               aes(xmin = xmin_bar, xmax = xmax_bar,
                   label = p.signif, y.position = y.position),
               inherit.aes = FALSE,
               tip.length  = 0.02, bracket.nudge.y = 0,
               size = 1.3, color = "grey20",
               label.size  = 14, fontface = "bold",
               vjust = -0.3) +

  scale_fill_gradient2(low = "#2980B9", mid = "#F0E442", high = "#C0392B",
                       midpoint = median(df_E$PTC.2.EJC), guide = "none") +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, Y_MAX_E), expand = c(0, 0)) +

  labs(x        = "PTC-to-EJC distance quintile (nt)",
       y        = "Mean NMD efficiency",
       title    = "NMD Efficiency by Distance Quintile",
       # plain ASCII: ggsave's default pdf() device cannot draw the "<=" glyph
       subtitle = if (E_MATCH_D) "Same penultimate exon variants as D"
                  else sprintf("All penultimate exon variants (n = %d)", nrow(df_E)),
       tag      = "E") +
  scatter_theme +
  theme(
    # 3-line labels (Q#, range, n), horizontal and centred under each bar.
    # With E matched to D the widest range is ~"200–250 nt", which fits.
    axis.text.x = element_text(size = FONT - 10, lineheight = 0.95,
                               face = "bold", angle = 0,
                               hjust = 0.5, vjust = 1)
  )

# PANELS F, G — penultimate exon length groups (same style as D)
pF <- make_scatter_panel(
  df_penu %>% filter(exon_len_grp == "0–200", PTC.2.EJC <= 200),
  col_medium, 200, seq(0, 200, 50),
  title = "Penultimate exon 0–200 bp", tag = "F", span = 0.80)

pG <- make_scatter_panel(
  df_penu %>% filter(exon_len_grp == ">200", PTC.2.EJC <= 400),
  col_long, 400, seq(0, 400, 100),
  title = "Penultimate exon >200 bp", tag = "G", span = 0.80)

# ══════════════════════════════════════════════════════════════════════════════
# COMBINE ALL PANELS
#   Row 1: A | B | C   (shared violin + box style; widths 3:2:2)
#   Row 2: D (scatter all penu, with inflection) | E (quintile bar all penu)
#   Row 3: F | G     (scatter per exon length group: 0-200 bp | >200 bp)
# ══════════════════════════════════════════════════════════════════════════════
row1 <- (pA | pB | pC_ejc) + plot_layout(widths = c(3, 2, 2))  # equal width per violin

combined <- row1 / (pD | pE) / (pF | pG) +
  plot_layout(heights = c(1, 1, 1)) &
  theme(plot.margin = margin(10, 14, 10, 14))

# Page width reduced (was 48) since the F row now has 2 panels instead of 4.
ggsave("Figure3_revised.pdf", combined,
       width = 36, height = 44, dpi = 300,
       limitsize = FALSE)
ggsave("Figure3_revised.png", combined,
       width = 36, height = 44, dpi = 300,
       limitsize = FALSE)

print(combined)
cat("\nSaved: Figure3_revised.pdf / .png\n")
