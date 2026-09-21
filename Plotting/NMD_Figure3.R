# ══════════════════════════════════════════════════════════════════════════════
# NMD — Merged Figure: Figure 3 A+B + Penultimate Exon PTC-to-EJC Analysis
# (FILTERED VERSION v5 — restricted to variants with ALLELE.RAT >= 0.35,
#                       penultimate set further restricted to transcripts with
#                       >= 3 exons, exon-length groups COLLAPSED to 2 tiers
#                       (0-200 bp vs >200 bp), per-group quintile bars (G row)
#                       REMOVED, global font size and line weights increased,
#                       inflection point on Panel D made more prominent)
#
# Layout:
#   Row 1: A (Exon position raincloud) | B (EJC downstream violin)
#                                       | C (EJC overlap violin)
#   Row 2: D (Scatter + loess + inflection point, all penultimate, x<=250)
#           | E (Quintile bar, all penultimate)
#   Row 3: F1 (0-200 bp scatter) | F2 (>200 bp scatter)
#
# Required file:
#   TOPMed_stopgain_withfullpenultimateMarch30.csv
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
library(ggdist)
library(gghalves)
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
df_raw <- read.csv("TOPMed_stopgain_withfullpenultimateMarch30.csv",
                   stringsAsFactors = FALSE) %>%
  filter(!is.na(ALLELE.RAT),
         ALLELE.RAT >= 0.35) %>%   # <<< NEW: restrict to variants with allele.rat >= 0.35
  mutate(
    # Panel A — exon position
    exon_pos = case_when(
      last.exon == "lastexon" ~ "Last exon",
      last.exon != "lastexon" &
        penultimate.exon == "penultimate.last50bp" ~ "Penultimate",
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
# PANEL A — Raincloud: PTC-bearing exon position (variant-level)
# ══════════════════════════════════════════════════════════════════════════════
pal_A <- c("Last exon"   = col_blue,
           "Penultimate" = col_grey,
           "Upstream"    = col_teal)

nA <- df_raw %>% count(exon_pos) %>%
  mutate(label = paste0("n = ", format(n, big.mark = ",")))

stat_A <- df_raw %>%
  wilcox_test(ALLELE.RAT ~ exon_pos,
              comparisons = list(c("Last exon","Penultimate"),
                                 c("Penultimate","Upstream"),
                                 c("Last exon","Upstream"))) %>%
  add_significance("p") %>%
  add_xy_position(x = "exon_pos") %>%
  mutate(y.position = c(1.08, 1.18, 1.28))

pA <- ggplot(df_raw, aes(x = exon_pos, y = ALLELE.RAT,
                          fill = exon_pos, color = exon_pos)) +

  stat_halfeye(adjust = 0.60, width = 0.50, justification = -0.28,
               .width = 0, point_colour = NA, alpha = 0.65) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 1.6,
               fill = "white", alpha = 0.88, color = "grey15") +
  geom_half_point(side = "l", range_scale = 0.28,
                  alpha = 0.14, size = 0.85, shape = 16) +
  stat_summary(fun = median, geom = "point", size = 4.5,
               shape = 21, fill = "white", color = "black", stroke = 1.3) +
  ref05 +
  geom_text(data = nA, aes(x = exon_pos, y = 1.42, label = label),
            inherit.aes = FALSE, size = 14,
            fontface = "bold", color = "grey15") +
  stat_pvalue_manual(stat_A, label = "p.signif",
                     tip.length = 0.012, bracket.size = 1.4,
                     size = 18, color = "grey20") +
  scale_fill_manual(values = pal_A) +
  scale_color_manual(values = pal_A) +
  y_nmd_dist +
  labs(y = "NMD efficiency", tag = "A", title = "Exon position") +
  base_theme +
  theme(
    # "Last exon" / "Penultimate" / "Upstream" are wide labels relative to
    # the column width — shrink + tilt slightly so they don't collide.
    axis.text.x = element_text(size = FONT - 4, face = "bold",
                               color = "grey10",
                               angle = 18, hjust = 1, vjust = 1)
  )

# ══════════════════════════════════════════════════════════════════════════════
# PANEL B — Violin + box: EJC downstream No/Yes (variant-level)
# ══════════════════════════════════════════════════════════════════════════════
nB <- df_raw %>% drop_na(ejc) %>% count(ejc) %>%
  mutate(label = paste0("n = ", format(n, big.mark = ",")))

stat_B <- df_raw %>% drop_na(ejc) %>%
  wilcox_test(ALLELE.RAT ~ ejc) %>%
  add_significance("p") %>% add_xy_position(x = "ejc") %>%
  mutate(y.position = 1.18)

pB <- ggplot(df_raw %>% drop_na(ejc),
             aes(x = ejc, y = ALLELE.RAT, fill = ejc, color = ejc)) +

  geom_violin(alpha = 0.42, width = 0.78, trim = TRUE, linewidth = 1.2) +
  geom_boxplot(width = 0.20, outlier.shape = 21, outlier.size = 2.5,
               outlier.alpha = 0.25, linewidth = 1.6,
               fill = "white", alpha = 0.88, color = "grey10") +
  stat_summary(fun = median, geom = "point", size = 4.5,
               shape = 21, fill = "white", color = "black", stroke = 1.3) +
  ref05 +
  geom_text(data = nB, aes(x = ejc, y = 1.42, label = label),
            inherit.aes = FALSE, size = 14,
            fontface = "bold", color = "grey15") +
  stat_pvalue_manual(stat_B, label = "p.signif",
                     tip.length = 0.012, bracket.size = 1.4,
                     size = 18, color = "grey20") +
  scale_fill_manual(values = pal_ejc) +
  scale_color_manual(values = pal_ejc) +
  y_nmd_dist +
  labs(y = "NMD efficiency", tag = "B", title = "EJC downstream of PTC") +
  base_theme

# ══════════════════════════════════════════════════════════════════════════════
# PANEL C — Boxplot: EJC Overlap (No vs Yes) — matches reference image
# ══════════════════════════════════════════════════════════════════════════════
df_raw <- df_raw %>%
  mutate(
    ejc_overlap = case_when(
      has_ejc_overlap == "yes" ~ "Yes",
      has_ejc_overlap == "no"  ~ "No",
      TRUE ~ NA_character_
    ),
    ejc_overlap = factor(ejc_overlap, levels = c("No", "Yes"))
  )

dfC_ejc <- df_raw %>% drop_na(ejc_overlap)

nC_ejc <- dfC_ejc %>% count(ejc_overlap) %>%
  mutate(label = format(n, big.mark = ","))

stat_C_ejc <- dfC_ejc %>%
  wilcox_test(ALLELE.RAT ~ ejc_overlap) %>%
  add_significance("p") %>% add_xy_position(x = "ejc_overlap") %>%
  mutate(y.position = 1.18)

pal_ejc_ov <- c("No" = "#2ECC71", "Yes" = "#9B59B6")

pC_ejc <- ggplot(dfC_ejc,
                  aes(x = ejc_overlap, y = ALLELE.RAT,
                      fill = ejc_overlap, color = ejc_overlap)) +

  geom_violin(alpha = 0.42, width = 0.78, trim = TRUE, linewidth = 1.2) +
  geom_boxplot(width = 0.22, outlier.shape = 21, outlier.size = 2.5,
               outlier.alpha = 0.28, linewidth = 1.6,
               fill = "white", alpha = 0.88, color = "grey10") +
  stat_summary(fun = median, geom = "point", size = 5,
               shape = 21, fill = "white", color = "black", stroke = 1.4) +
  ref05 +
  geom_text(data = nC_ejc,
            aes(x = ejc_overlap, y = 1.42, label = label),
            inherit.aes = FALSE, size = 14,
            fontface = "bold", color = "grey15") +
  stat_pvalue_manual(stat_C_ejc, label = "p.signif",
                     tip.length = 0.012, bracket.size = 1.4,
                     size = 18, color = "grey20") +
  scale_fill_manual(values  = pal_ejc_ov) +
  scale_color_manual(values = pal_ejc_ov) +
  y_nmd_dist +
  labs(y = "NMD efficiency", tag = "C",
       title = "EJC Overlap") +
  base_theme

# ══════════════════════════════════════════════════════════════════════════════
# PANEL D — Scatter + loess + INFLECTION POINT (all penultimate, x <= 250)
# ══════════════════════════════════════════════════════════════════════════════
X_MAX_D  <- 250
n_rm_D   <- sum(df_penu$PTC.2.EJC > X_MAX_D)
df_D     <- df_penu %>% filter(PTC.2.EJC <= X_MAX_D)
r_D      <- cor(df_D$PTC.2.EJC, df_D$ALLELE.RAT, use = "complete.obs",
                method = "spearman")   # <<< v5: Spearman (was default Pearson)
minima_D <- find_local_minima(df_D$PTC.2.EJC, df_D$ALLELE.RAT, span = 0.75)

cat(sprintf("\nPanel D: n=%d (removed %d >%d nt), r=%.3f\n",
            nrow(df_D), n_rm_D, X_MAX_D, r_D))

pD <- ggplot(df_D, aes(x = PTC.2.EJC, y = ALLELE.RAT)) +

  geom_point(aes(color = PTC.2.EJC), size = 7.5, alpha = 0.65, shape = 16) +
  geom_smooth(method = "loess", span = 0.75, color = "black",
              linewidth = 3.2, fill = "grey70", alpha = 0.25, se = TRUE) +

  {if (nrow(minima_D) > 0) list(
    geom_segment(data = minima_D,
                 aes(x = x_min, xend = x_min, y = 0, yend = y_min),
                 color = MIN_COL, linetype = "dashed",
                 linewidth = 1.7, inherit.aes = FALSE),
    geom_point(data = minima_D,
               aes(x = x_min, y = y_min),
               shape = 25, size = 11, fill = MIN_COL, color = "white",
               stroke = 1.6, inherit.aes = FALSE),
    geom_label(data = minima_D,
               aes(x = x_min + 4, y = y_min - 0.13,
                   label = sprintf("Inflection\n(%.0f nt)", x_min)),
               hjust = 0, size = 16, fontface = "bold.italic",
               color = MIN_COL, fill = alpha("white", 0.92),
               label.size = 0.8, inherit.aes = FALSE)
  )} +

  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "grey40", linewidth = 1.04) +
  annotate("text", x = 155, y = 0.08,
           label = sprintf("r = %.3f", r_D),
           size = 18, fontface = "bold.italic", color = "grey20") +
  annotate("text", x = 3, y = 0.98,
           label = paste0("n = ", nrow(df_D)),
           hjust = 0, size = 18, fontface = "bold", color = "grey20") +

  scale_color_gradient2(low = "#2980B9", mid = "#F0E442", high = "#C0392B",
                        midpoint = 125, name = "PTC-to-EJC\n(nt)",
                        guide = guide_colorbar(barwidth = 1.8, barheight = 14,
                                               title.hjust = 0.5)) +
  scale_x_continuous(breaks = seq(0, X_MAX_D, 50),
                     limits = c(-2, X_MAX_D + 8), expand = c(0.01, 0)) +
  y_nmd_scat +

  labs(x        = "PTC-to-EJC distance (nt)",
       y        = "NMD efficiency",
       title    = "PTC-to-EJC Distance vs NMD Efficiency",
       subtitle = "All penultimate exon variants",
       tag      = "D") +
  scatter_theme +
  theme(legend.position = "right",
        legend.title    = element_text(size = FONT - 1, face = "bold"),
        legend.text     = element_text(size = FONT - 2))

# ══════════════════════════════════════════════════════════════════════════════
# PANEL E — Quintile bar (all penultimate variants)
# ══════════════════════════════════════════════════════════════════════════════
bin_E <- df_penu %>%
  mutate(
    dist_q = cut(PTC.2.EJC,
                  breaks = quantile(PTC.2.EJC, probs = seq(0, 1, 0.2),
                                    na.rm = TRUE),
                  include.lowest = TRUE, labels = FALSE)
  ) %>%
  drop_na(dist_q) %>%
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
  mutate(lab = sprintf("Q%d\n%d\u2013%d nt\n(n=%d)",
                        dist_q, round(d_min), round(d_max), n))

# Sequential Wilcoxon tests between adjacent quintiles
stat_E_bar <- df_penu %>%
  mutate(dist_q_fac = factor(cut(PTC.2.EJC,
                                  breaks = quantile(PTC.2.EJC, probs = seq(0, 1, 0.2),
                                                    na.rm = TRUE),
                                  include.lowest = TRUE, labels = FALSE))) %>%
  drop_na(dist_q_fac) %>%
  wilcox_test(ALLELE.RAT ~ dist_q_fac,
              comparisons = list(c("1","2"),c("2","3"),c("3","4"),c("4","5"))) %>%
  add_significance("p") %>%
  add_xy_position(x = "dist_q_fac") %>%
  mutate(y.position = c(0.82, 0.87, 0.82, 0.87))

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
  geom_text(aes(label = sprintf("%.3f", mean_NMD), y = ci_hi + 0.018),
            size = 14, fontface = "bold", color = "grey20") +

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
                       midpoint = median(df_penu$PTC.2.EJC), guide = "none") +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 0.98), expand = c(0, 0)) +

  labs(x        = "PTC-to-EJC distance quintile",
       y        = "Mean NMD efficiency",
       title    = "NMD Efficiency by Distance Quintile",
       subtitle = "All penultimate exon variants",
       tag      = "E") +
  scatter_theme +
  theme(
    # Each bar's label is 3 lines (Q#, range, n=…). The widest range
    # ("129–4291 nt") overflows the column width — shrink the font
    # and tilt slightly so the labels stop colliding.
    axis.text.x = element_text(size = FONT - 4, lineheight = 1.0,
                               face = "bold", angle = 12,
                               hjust = 1, vjust = 1)
  )

# ══════════════════════════════════════════════════════════════════════════════
# PANELS F1/F2 — scatter + loess per exon length group (2 groups)
#   0-200 bp        | xmax = 200
#   >200 bp         | xmax = 1000
# ══════════════════════════════════════════════════════════════════════════════
grp_cfg <- list(
  `0–200` = list(col = col_medium, xmax = 200, xbrk = seq(0, 200, 50)),
  `>200`  = list(col = col_long,   xmax = 400, xbrk = seq(0, 400, 100))
)
grp_titles_F2 <- c(
  `0–200` = "Penultimate exon 0–200 bp",
  `>200`  = "Penultimate exon >200 bp"
)

make_F2_panel <- function(grp) {
  cfg <- grp_cfg[[grp]]
  sub <- df_penu %>% filter(exon_len_grp == grp, PTC.2.EJC <= cfg$xmax)
  if (nrow(sub) < 5) {
    # Too few points to fit loess / compute correlation reliably
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = sprintf("%s\n(n = %d, too few)",
                                 grp_titles_F2[grp], nrow(sub)),
                 size = 8, color = "grey30") +
        theme_void()
    )
  }
  r_g <- cor(sub$PTC.2.EJC, sub$ALLELE.RAT, use = "complete.obs",
             method = "spearman")   # <<< v5: Spearman (was default Pearson)

  ggplot(sub, aes(x = PTC.2.EJC, y = ALLELE.RAT)) +
    geom_point(color = cfg$col, size = 7.5, alpha = 0.65, shape = 16) +
    geom_smooth(method = "loess", span = 0.80, color = "black",
                linewidth = 3.0, fill = "grey70", alpha = 0.22, se = TRUE) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "grey40", linewidth = 1.04) +
    annotate("text", x = cfg$xmax * 0.58, y = 0.08,
             label = sprintf("r = %.3f", r_g),
             size = 18, fontface = "bold.italic", color = "grey20") +
    annotate("text", x = cfg$xmax * 0.02 + 1, y = 0.98,
             label = paste0("n = ", nrow(sub)),
             hjust = 0, size = 18, fontface = "bold", color = "grey20") +
    scale_x_continuous(breaks = cfg$xbrk,
                       limits = c(-1, cfg$xmax + cfg$xmax * 0.04),
                       expand = c(0.01, 0)) +
    y_nmd_scat +
    labs(x = "PTC-to-EJC distance (nt)", y = "NMD efficiency",
         title = grp_titles_F2[grp]) +
    scatter_theme +
    theme(plot.title   = element_text(size = FONT + 2, face = "bold",
                                      color = TITLE_COL, hjust = 0.5),
          plot.caption = element_blank(),
          plot.margin  = margin(10, 14, 10, 14))
}

pF1 <- make_F2_panel("0\u2013200")
pF2 <- make_F2_panel(">200")

pF_row <- (pF1 | pF2) +
  plot_annotation(
    tag_levels = list(c("F","")),
    theme = theme(
      plot.tag     = element_text(size = FONT + 28, face = "bold")
    )
  )

# ══════════════════════════════════════════════════════════════════════════════
# COMBINE ALL PANELS
#   Row 1: A (raincloud exon position) | B (violin EJC) | C (violin EJC overlap)
#   Row 2: D (scatter all penu, with inflection) | E (quintile bar all penu)
#   Row 3: F1 | F2   (scatter per exon length group: 0-200 bp | >200 bp)
# ══════════════════════════════════════════════════════════════════════════════
combined <- (pA | pB | pC_ejc) / (pD | pE) / pF_row +
  plot_layout(heights = c(1, 1, 1)) &
  theme(plot.margin = margin(10, 14, 10, 14))

# Page width reduced (was 48) since the F row now has 2 panels instead of 4.
ggsave("NMD_Figure_merged_filtered_AR0.35.pdf", combined,
       width = 36, height = 44, dpi = 300,
       limitsize = FALSE)
ggsave("NMD_Figure_merged_filtered_AR0.35.png", combined,
       width = 36, height = 44, dpi = 300,
       limitsize = FALSE)

print(combined)
cat("\nSaved: NMD_Figure_merged_filtered_AR0.35.pdf / .png\n")
