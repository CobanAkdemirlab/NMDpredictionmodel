# ══════════════════════════════════════════════════════════════════════════════
# Supplemental Figure 6 — relative PTC location x transcript architecture
#
# Panels (2 x 2):
#   A = NMD efficiency by relative-PTC-location tertile, faceted by CDS length
#   B = Relative PTC location vs NMD efficiency, LOESS per CDS-length group
#   C = As A, faceted by exon count
#   D = As B, one LOESS per exon-count group
#
#   Points are medians of ALLELE.RAT and relativePTClocation over a unit of
#   variants. The unit and the Low/Medium/High cut-offs are chosen to
#   reproduce the original figure (see "RECOVER THE ORIGINAL METHOD"); the
#   choice is printed and should be stated in the legend.
#   Tests: pairwise Wilcoxon rank-sum (Low-Medium, Medium-High, Low-High)
#   within each facet, Benjamini-Hochberg adjusted.
#
# Input : TOPMed_stopgain_September25_corrected_readyformodel.csv
# Filter: ALLELE.RAT >= 0.35  (same as Figures 3, 4, 6 and Supplemental Figure 2)
# Output: SupplementalFigure4_revised.pdf / .png
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(patchwork)
library(rstatix)
library(ggpubr)

set.seed(42)

# ── Config ───────────────────────────────────────────────────────────────────
INPUT_CSV <- "TOPMed_stopgain_September25_corrected_readyformodel.csv"
AR_MIN    <- 0.35
CDS_CUTS  <- c(1289, 2181)
EXON_CUTS <- c(5, 10, 20)

TITLE_COL <- "#2E6DA4"
FONT      <- 22
PAL_TERT  <- c("Low" = "#7FBF7B", "Medium" = "#6A9FD8", "High" = "#E07B7B")
PAL_CDS   <- c("Short (<1,289 bp)"            = "#E8924A",
               "Medium (1,289–2,181 bp)" = "#6A9FD8",
               "Long (>2,181 bp)"             = "#3EC9A7")
PAL_EXON  <- c("2–5 exons"   = "#D9534F", "6–10 exons" = "#E8924A",
               "11–20 exons" = "#4A72B0", ">20 exons"       = "#3EC9A7")

# ── Load & filter ───────────────────────────────────────────────────────────
df <- read.csv(INPUT_CSV, stringsAsFactors = FALSE) %>%
  filter(!is.na(ALLELE.RAT), ALLELE.RAT >= AR_MIN)
cat(sprintf("[Filter: ALLELE.RAT >= %.2f] variants: %d (main figures: 5,749)\n",
            AR_MIN, nrow(df)))

required <- c("TxName", "relativePTClocation", "cds_length", "exon_count")
missing  <- setdiff(required, names(df))
if (length(missing) > 0) stop("Missing column(s): ", paste(missing, collapse = ", "))

df <- df %>%
  drop_na(relativePTClocation) %>%
  mutate(
    cds_grp = case_when(cds_length <  CDS_CUTS[1] ~ "Short (<1,289 bp)",
                        cds_length <= CDS_CUTS[2] ~ "Medium (1,289–2,181 bp)",
                        TRUE                      ~ "Long (>2,181 bp)"),
    cds_grp = factor(cds_grp, levels = names(PAL_CDS)),
    exon_grp = case_when(exon_count <= EXON_CUTS[1] ~ "2–5 exons",
                         exon_count <= EXON_CUTS[2] ~ "6–10 exons",
                         exon_count <= EXON_CUTS[3] ~ "11–20 exons",
                         TRUE                       ~ ">20 exons"),
    exon_grp = factor(exon_grp, levels = names(PAL_EXON)))
if ("downstream" %in% names(df))
  df <- df %>% mutate(ejc = case_when(downstream == "upstream of last EJC"   ~ "Yes",
                                      downstream == "downstream of last EJC" ~ "No"))

# ══════════════════════════════════════════════════════════════════════════════
# RECOVER THE ORIGINAL METHOD
#   The original (June) S4 had ~4,528 points in panel A (1,657 / 1,500 / 1,371
#   per CDS group), i.e. one point per transcript AND location group: the facet
#   totals of that unit match the original within 1-3%. The split into
#   Low / Medium / High, however, differed by facet (e.g. >20 exons: original
#   Low = 146 vs 203 with global quartiles, while 2-5 exons matched), which a
#   single global cut-off cannot produce. So the cut-offs may have been computed
#   WITHIN each facet (CDS-length group in A, exon-count group in C).
#
#   The script tries every combination of
#     unit : tx (one point per transcript), tx_cat (per transcript x location
#            group), tx_ejc_cat (per transcript x EJC status x location group)
#     rule : fixed cut-offs | global quantiles | quantiles within each facet,
#            as tertiles (1/3, 2/3) or quartile groups (25%, 75%)
#   builds panels A and C separately (their facets differ), compares the counts
#   with the original figure, and uses the closest combination (printed).
#   Set UNIT_KEY / RULE_KEY to force a choice.
# ══════════════════════════════════════════════════════════════════════════════
UNIT_KEY <- "tx_cat"
RULE_KEY <- "global_tertiles"

# Original counts: panel A (CDS Short/Medium/Long x Low/Medium/High),
# panel C (2-5 / 6-10 / 11-20 / >20 exons x Low/Medium/High)
OLD_A <- c(505, 735, 417,  384, 745, 371,  299, 715, 357)
OLD_C <- c(258, 417, 227,  381, 648, 344,  403, 737, 367,  146, 395, 207)

TERT  <- c(1/3, 2/3); QUART <- c(0.25, 0.75)
rules <- list(
  fixed_fig4         = list(type = "fixed",  v = c(0.3750, 0.7144)),
  fixed_25_75        = list(type = "fixed",  v = c(0.25, 0.75)),
  fixed_thirds       = list(type = "fixed",  v = TERT),
  global_tertiles    = list(type = "global", v = TERT),
  global_quartiles   = list(type = "global", v = QUART),
  facet_tertiles     = list(type = "facet",  v = TERT),
  facet_quartiles    = list(type = "facet",  v = QUART)
)
cut3 <- function(x, cuts)
  factor(case_when(x <= cuts[1] ~ "Low", x <= cuts[2] ~ "Medium", TRUE ~ "High"),
         levels = c("Low", "Medium", "High"))
# category for values x, with facet labels g (used only by "facet" rules)
cat_by <- function(x, g, rule) {
  if (rule$type == "fixed")  return(cut3(x, rule$v))
  if (rule$type == "global") return(cut3(x, quantile(x, rule$v, na.rm = TRUE)))
  out <- factor(rep(NA, length(x)), levels = c("Low", "Medium", "High"))
  for (k in unique(g[!is.na(g)])) {
    i <- which(g == k)
    out[i] <- cut3(x[i], quantile(x[i], rule$v, na.rm = TRUE))
  }
  out
}

unit_sets <- list(
  tx         = list(keys = "TxName",           split = FALSE),
  tx_cat     = list(keys = "TxName",           split = TRUE),
  tx_ejc_cat = list(keys = c("TxName", "ejc"), split = TRUE)
)
unit_sets <- Filter(function(u) all(u$keys %in% names(df)), unit_sets)

summ <- function(d) summarise(d, ALLELE.RAT = median(ALLELE.RAT),
                              relativePTClocation = median(relativePTClocation),
                              cds_grp = first(cds_grp), exon_grp = first(exon_grp),
                              n_var = n(), .groups = "drop")
# Build the points for one panel; 'strat' = the facet variable of that panel
build_units <- function(u, rule, strat) {
  if (u$split) {
    df %>% mutate(ptc_cat = cat_by(relativePTClocation, .data[[strat]], rule)) %>%
      group_by(across(all_of(c(u$keys, "ptc_cat")))) %>% summ()
  } else {
    df %>% group_by(across(all_of(u$keys))) %>% summ() %>%
      mutate(ptc_cat = cat_by(relativePTClocation, .data[[strat]], rule))
  }
}
cnt <- function(t, strat) t %>% count(.data[[strat]], ptc_cat) %>%
  complete(.data[[strat]], ptc_cat, fill = list(n = 0)) %>%
  arrange(.data[[strat]], ptc_cat) %>% pull(n)

fit <- expand_grid(unit = names(unit_sets), rule = names(rules)) %>%
  mutate(res = map2(unit, rule, function(un, ru) {
    a <- cnt(build_units(unit_sets[[un]], rules[[ru]], "cds_grp"),  "cds_grp")
    c <- cnt(build_units(unit_sets[[un]], rules[[ru]], "exon_grp"), "exon_grp")
    tibble(total_A = sum(a),
           mean_abs_diff_vs_original = round(mean(abs(c(a, c) - c(OLD_A, OLD_C))), 1),
           A_counts = paste(a, collapse = "/"),
           C_counts = paste(c, collapse = "/"))
  })) %>% unnest(res) %>% arrange(mean_abs_diff_vs_original)

cat("\n── Method check ──\n")
cat("Original A: 505/735/417/384/745/371/299/715/357\n")
cat("Original C: 258/417/227/381/648/344/403/737/367/146/395/207\n")
print(as.data.frame(fit), row.names = FALSE)

unit_pick <- if (UNIT_KEY == "auto") fit$unit[1] else UNIT_KEY
rule_pick <- if (RULE_KEY == "auto") fit$rule[1] else RULE_KEY
RULE <- rules[[rule_pick]]
cat(sprintf("\nUsing unit: %s | rule: %s\n", unit_pick, rule_pick))
txA <- build_units(unit_sets[[unit_pick]], RULE, "cds_grp")
txC <- build_units(unit_sets[[unit_pick]], RULE, "exon_grp")
cat(sprintf("Points plotted: A/B %d, C/D %d\n", nrow(txA), nrow(txC)))

# print the actual cut-off values used (per facet for facet rules)
show_cuts <- function(strat) {
  x <- if (unit_sets[[unit_pick]]$split) df else df %>% group_by(TxName) %>% summ()
  if (RULE$type == "facet") {
    print(x %>% group_by(.data[[strat]]) %>%
            summarise(cut_low  = round(quantile(relativePTClocation, RULE$v[1]), 4),
                      cut_high = round(quantile(relativePTClocation, RULE$v[2]), 4),
                      .groups = "drop"))
  } else if (RULE$type == "global") {
    print(round(quantile(x$relativePTClocation, RULE$v), 4))
  } else print(RULE$v)
}
cat("Cut-offs, panel A:\n"); show_cuts("cds_grp")
cat("Cut-offs, panel C:\n"); show_cuts("exon_grp")

UNIT_DESC <- c(tx         = "one point per transcript",
               tx_cat     = "one point per transcript and location group",
               tx_ejc_cat = "one point per transcript, EJC status and location group")[[unit_pick]]
grp_word <- if (identical(RULE$v, QUART)) "quartile groups: Q1 | Q2\u2013Q3 | Q4" else "tertiles: lower | middle | upper third of variants"
XLAB <- switch(RULE$type,
  fixed  = sprintf("Relative PTC location (Low \u2264%.2f, Medium %.2f\u2013%.2f, High >%.2f)",
                   RULE$v[1], RULE$v[1], RULE$v[2], RULE$v[2]),
  global = if (identical(RULE$v, QUART)) sprintf("Relative PTC location (%s)", grp_word) else {
    qq <- quantile((if (unit_sets[[unit_pick]]$split) df else df %>% group_by(TxName) %>% summ())$relativePTClocation, RULE$v)
    sprintf("Relative PTC location tertile (\u2264%.2f | %.2f\u2013%.2f | >%.2f)", qq[1], qq[1], qq[2], qq[2]) },
  facet  = sprintf("Relative PTC location (%s, within each group)", grp_word))

# ── Theme ────────────────────────────────────────────────────────────────────
base_theme <- theme_classic(base_size = FONT) +
  theme(plot.title       = element_text(size = FONT + 2, face = "bold",
                                        color = TITLE_COL, hjust = 0.5),
        axis.title       = element_text(face = "bold"),
        axis.text        = element_text(color = "grey10"),
        strip.text       = element_text(size = FONT - 1, face = "bold",
                                        lineheight = 0.95),
        strip.background = element_rect(fill = "grey94", color = "grey60"),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
        plot.tag         = element_text(size = FONT + 12, face = "bold"),
        plot.margin      = margin(10, 14, 10, 14))

# ── Violin panel (A, C) ──────────────────────────────────────────────────────
make_violin <- function(d, strat, title, tag) {
  d <- d %>% mutate(strat = {{ strat }})
  n_lab <- d %>% count(strat, ptc_cat) %>% mutate(label = prettyNum(n, big.mark = ","))

  st <- d %>% group_by(strat) %>%
    wilcox_test(ALLELE.RAT ~ ptc_cat,
                comparisons = list(c("Low", "Medium"), c("Medium", "High"),
                                   c("Low", "High")),
                p.adjust.method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "ptc_cat") %>%
    mutate(y.position = rep(c(1.04, 1.12, 1.21), length.out = n()),
           ptc_cat = group1)   # needed by inherited aes

  cat(sprintf("\n[%s] %s\n", tag, title))
  print(n_lab %>% mutate(strat = gsub("\n", " ", strat)) %>% select(strat, ptc_cat, n) %>%
          pivot_wider(names_from = ptc_cat, values_from = n))
  print(st %>% mutate(strat = gsub("\n", " ", strat)) %>%
          select(strat, group1, group2, p, p.adj, p.adj.signif))

  ggplot(d, aes(x = ptc_cat, y = ALLELE.RAT, fill = ptc_cat, color = ptc_cat)) +
    geom_violin(alpha = 0.35, trim = TRUE, linewidth = 0.6, width = 0.9) +
    geom_jitter(width = 0.15, height = 0, size = 0.5, alpha = 0.25, shape = 16) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
                 color = "grey10", linewidth = 0.6) +
    stat_summary(fun = median, geom = "point", shape = 23, size = 2.8,
                 fill = "white", color = "grey10") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50",
               linewidth = 0.6) +
    geom_text(data = n_lab, aes(x = ptc_cat, y = 1.37, label = label),
              inherit.aes = FALSE, size = 5.5, fontface = "bold", color = "grey10") +
    stat_pvalue_manual(st, label = "p.adj.signif", tip.length = 0.01,
                       bracket.size = 0.6, size = 6, color = "grey20") +
    facet_wrap(~ strat, nrow = 1,
               labeller = as_labeller(function(x)            # two-line strips
                 sub(" \\(", "\n(", sub(" exons", "\nexons", x)))) +
    scale_fill_manual(values = PAL_TERT) +
    scale_color_manual(values = PAL_TERT) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(-0.02, 1.42),
                       expand = c(0, 0)) +
    labs(x = XLAB, y = "NMD efficiency",
         title = title, tag = tag) +
    base_theme +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))  # 4 facets in C
}

# ── LOESS panel (B, D) ───────────────────────────────────────────────────────
make_loess <- function(d, strat, pal, legend_title, title, tag) {
  ggplot(d, aes(x = relativePTClocation, y = ALLELE.RAT,
                color = {{ strat }}, fill = {{ strat }})) +
    geom_point(size = 1.2, alpha = 0.35, shape = 16) +
    geom_smooth(method = "loess", formula = y ~ x, span = 0.75,
                linewidth = 1.6, alpha = 0.15, se = TRUE) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50",
               linewidth = 0.6) +
    scale_color_manual(values = pal, name = legend_title) +
    scale_fill_manual(values = pal, name = legend_title) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1),
                       expand = c(0.01, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.02),
                       expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = "Relative PTC location", y = "NMD efficiency",
         title = title, tag = tag) +
    base_theme +
    theme(legend.position = "right",
          legend.title = element_text(face = "bold"),
          legend.text  = element_text(size = FONT - 4),
          legend.key.height = unit(1.4, "lines"))
}

pA <- make_violin(txA, cds_grp,  "NMD efficiency by CDS length", "A")
pB <- make_loess(txA, cds_grp, PAL_CDS, "CDS length",
                 "PTC location vs NMD by CDS length", "B")
pC <- make_violin(txC, exon_grp, "NMD efficiency by exon count", "C")
pD <- make_loess(txC, exon_grp, PAL_EXON, "Exon count",
                 "PTC location vs NMD by exon count", "D")

final <- (pA | pB) / (pC | pD) +
  plot_layout(widths = c(1.35, 1)) +
  plot_annotation(
    caption = sprintf(paste0("Filter: ALLELE.RAT >= %.2f  |  %s (median ALLELE.RAT and ",
                             "median relativePTClocation)  |  Wilcoxon tests, BH-adjusted"),
                      AR_MIN, UNIT_DESC),
    theme = theme(plot.caption = element_text(size = FONT - 6, color = "grey45",
                                              face = "italic", hjust = 0.5)))

ggsave("SupplementalFigure4_revised.pdf", final, width = 26, height = 15,
       limitsize = FALSE)
ggsave("SupplementalFigure4_revised.png", final, width = 26, height = 15,
       dpi = 300, limitsize = FALSE)
cat("\nSaved: SupplementalFigure4_revised.pdf / .png\n")
