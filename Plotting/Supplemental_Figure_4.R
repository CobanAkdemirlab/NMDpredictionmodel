# ══════════════════════════════════════════════════════════════════════════════
# NMD Supplemental Figure 4 — downstream-EJC effect across strata
#   (rebuilt October 2026 for the corrected dataset; original script lost)
#
# Panels (2 x 2):
#   A = Gene tolerance x EJC    (Highly intolerant | Highly tolerant)  transcript-level
#   B = Allele frequency x EJC  (Ultra-rare | Rare/Common)             variant-level
#   C = CDS length x EJC        (Short | Medium | Long)                transcript-level
#   D = Exon count x EJC        (2-5 | 6-10 | 11-20 | >20)             transcript-level
#
#   "EJC" = EJC downstream of the PTC (No / Yes), defined exactly as in
#   Figure 3B: downstream == "upstream of last EJC" -> Yes,
#              downstream == "downstream of last EJC" -> No.
#
#   Transcript-level = median ALLELE.RAT per transcript within each EJC group
#   (one point per TxName x EJC status). Allele frequency is a property of the
#   variant, not the transcript, so panel B stays variant-level.
#
#   Test: Wilcoxon rank-sum, No vs Yes within each stratum.
#
# Input : TOPMed_stopgain_September25_corrected_readyformodel.csv
# Filter: ALLELE.RAT >= 0.35  (same as Figures 3, 4, 6 and Supplemental Figure 2)
# Output: SupplementalFigure3_revised.pdf / .png
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(patchwork)
library(rstatix)
library(ggpubr)

set.seed(42)

# ── Config ───────────────────────────────────────────────────────────────────
INPUT_CSV <- "TOPMed_stopgain_September25_corrected_readyformodel.csv"
AR_MIN    <- 0.35
CDS_CUTS  <- c(1289, 2181)          # bp, as in Figure 4 / Supplemental Figure 2
EXON_CUTS <- c(5, 10, 20)

TITLE_COL <- "#2E6DA4"
FONT      <- 22
PAL_EJC   <- c("No" = "#E8924A", "Yes" = "#3EC9A7")   # Figure 3B colours

# ── Load & filter ────────────────────────────────────────────────────────────
df <- read.csv(INPUT_CSV, stringsAsFactors = FALSE) %>%
  filter(!is.na(ALLELE.RAT), ALLELE.RAT >= AR_MIN)
cat(sprintf("[Filter: ALLELE.RAT >= %.2f] variants: %d (main figures: 5,749)\n",
            AR_MIN, nrow(df)))

# Gene tolerance: use pLI.cat; rebuild from numeric pLI if absent
if (!"pLI.cat" %in% names(df)) {
  pli_col <- intersect(c("pLI", "pli", "gnomad_pLI", "pLI_score"), names(df))
  if (length(pli_col) == 0) stop("Neither 'pLI.cat' nor a numeric pLI column found.")
  cat(sprintf("pLI.cat not found; rebuilding from '%s'\n", pli_col[1]))
  df$pLI.cat <- case_when(
    df[[pli_col[1]]] >= 0.65 ~ "highly intolerant (pLI>=0.65)",
    df[[pli_col[1]]] <  0.35 ~ "highly tolerant (pLI<0.35)",
    TRUE ~ NA_character_)
}
required <- c("TxName", "downstream", "pLI.cat", "Freq.cat", "cds_length", "exon_count")
missing  <- setdiff(required, names(df))
if (length(missing) > 0) stop("Missing column(s): ", paste(missing, collapse = ", "))

df <- df %>%
  mutate(
    ejc = case_when(downstream == "upstream of last EJC"   ~ "Yes",
                    downstream == "downstream of last EJC" ~ "No"),
    ejc = factor(ejc, levels = c("No", "Yes")),
    gene_tol = case_when(
      pLI.cat == "highly intolerant (pLI>=0.65)" ~ "Highly\nintolerant",
      pLI.cat == "highly tolerant (pLI<0.35)"    ~ "Highly\ntolerant"),
    gene_tol = factor(gene_tol, levels = c("Highly\nintolerant", "Highly\ntolerant")),
    af_grp = case_when(Freq.cat == "Ultra-rare variants"  ~ "Ultra-rare",
                       Freq.cat == "Rare/Common variants" ~ "Rare/Common"),
    af_grp = factor(af_grp, levels = c("Ultra-rare", "Rare/Common")),
    cds_grp = case_when(
      cds_length <  CDS_CUTS[1] ~ "Short\n(<1,289 bp)",
      cds_length <= CDS_CUTS[2] ~ "Medium\n(1,289–2,181 bp)",
      TRUE                      ~ "Long\n(>2,181 bp)"),
    cds_grp = factor(cds_grp, levels = c("Short\n(<1,289 bp)",
                                         "Medium\n(1,289–2,181 bp)",
                                         "Long\n(>2,181 bp)")),
    exon_grp = case_when(
      exon_count <= EXON_CUTS[1] ~ "2–5\nexons",
      exon_count <= EXON_CUTS[2] ~ "6–10\nexons",
      exon_count <= EXON_CUTS[3] ~ "11–20\nexons",
      TRUE                       ~ ">20\nexons"),
    exon_grp = factor(exon_grp, levels = c("2–5\nexons", "6–10\nexons",
                                           "11–20\nexons", ">20\nexons"))
  )

cat("EJC groups:\n"); print(table(df$ejc, useNA = "ifany"))

# ── Transcript-level summary: median ALLELE.RAT per transcript x EJC group ──
# (stratum variables are transcript properties, so first() is safe)
#
# WHICH UNIT DID THE ORIGINAL FIGURE USE?
#   The original S3 had 4,295 points in panels C and D (No 1,200 / Yes 3,095),
#   but the original S2 had only 3,212 transcripts in total, so the original
#   "transcript-level" points were NOT one per transcript x EJC (that gives
#   ~3,700). The original code most likely built ONE summary table grouped by
#   transcript, EJC status and every stratum variable at once; allele-frequency
#   class varies between variants of the same transcript, so that splits a
#   transcript into extra points.
#   Rather than guess, the script below collapses the data with each candidate
#   grouping, compares the panel C and D counts with the original figure, and
#   uses the closest one (printed). Set UNIT_KEY to a candidate name to force it.
UNIT_KEY <- "auto"

# Original (June) Figure S3 counts, panels C and D, No/Yes per stratum
OLD_C <- c(576, 1072, 392, 1044, 232, 979)
OLD_D <- c(484, 409, 359, 1009, 272, 1114, 85, 563)

# Candidate groupings (each also includes EJC status). Columns that are not in
# the file are skipped automatically.
cand <- list(
  tx_ejc            = c("TxName"),
  tx_ejc_af         = c("TxName", "af_grp"),
  tx_ejc_allstrata  = c("TxName", "gene_tol", "af_grp", "cds_grp", "exon_grp"),
  tx_ejc_exonpos    = c("TxName", "last.exon", "penultimate.last50bp"),
  tx_ejc_mutexon    = c("TxName", "mut.exon"),
  gene_ejc          = c("gene"),
  gene_ejc_af       = c("gene", "af_grp"),
  variant_site      = c("TxName", "PTC.2.start")   # one point per PTC position
)
cand <- Filter(function(k) all(k %in% names(df)), cand)

collapse_by <- function(data, keys) {
  data %>%
    group_by(across(all_of(c(keys, "ejc")))) %>%
    summarise(ALLELE.RAT = median(ALLELE.RAT),
              cds_grp    = first(cds_grp),  exon_grp = first(exon_grp),
              gene_tol   = first(gene_tol), af_grp   = first(af_grp),
              n_var      = n(), .groups = "drop")
}
counts_vec <- function(u, strat) {
  u %>% drop_na(ejc, {{ strat }}) %>% count({{ strat }}, ejc) %>%
    complete({{ strat }}, ejc, fill = list(n = 0)) %>% arrange({{ strat }}, ejc) %>% pull(n)
}

fit <- map_dfr(names(cand), function(nm) {
  u  <- collapse_by(df %>% drop_na(ejc), cand[[nm]])
  cC <- counts_vec(u, cds_grp); cD <- counts_vec(u, exon_grp)
  tibble(candidate = nm, keys = paste(cand[[nm]], collapse = " + "),
         total_C = sum(cC),
         mean_abs_diff_vs_original = mean(abs(c(cC, cD) - c(OLD_C, OLD_D))),
         C_counts = paste(cC, collapse = "/"))
}) %>% arrange(mean_abs_diff_vs_original)

cat("\n── Unit-of-analysis check (original S3: total 4,295; C = 576/1072/392/1044/232/979) ──\n")
print(as.data.frame(fit), row.names = FALSE)

chosen <- if (UNIT_KEY == "auto") fit$candidate[1] else UNIT_KEY
if (!chosen %in% names(cand)) stop("UNIT_KEY '", chosen, "' not available. Options: ",
                                   paste(names(cand), collapse = ", "))
cat(sprintf("\nUsing unit: %s  (%s + ejc)\n", chosen, paste(cand[[chosen]], collapse = " + ")))
UNITS <- collapse_by(df %>% drop_na(ejc), cand[[chosen]])

to_tx <- function(data, strat) {          # 'data' kept for call compatibility
  UNITS %>% drop_na(ejc, {{ strat }}) %>% mutate(strat = {{ strat }})
}
to_var <- function(data, strat) {
  data %>% drop_na(ejc, {{ strat }}) %>% mutate(strat = {{ strat }})
}

# ── Theme ────────────────────────────────────────────────────────────────────
base_theme <- theme_classic(base_size = FONT) +
  theme(legend.position  = "none",
        plot.title       = element_text(size = FONT + 2, face = "bold",
                                        color = TITLE_COL, hjust = 0.5),
        plot.caption     = element_text(size = FONT - 6, color = "grey45",
                                        hjust = 0.5),
        axis.title.x     = element_blank(),
        axis.title.y     = element_text(face = "bold"),
        axis.text.x      = element_text(size = FONT, face = "bold", color = "grey10"),
        strip.text       = element_text(size = FONT - 1, face = "bold",
                                        lineheight = 0.95),
        strip.background = element_rect(fill = "grey94", color = "grey60"),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
        plot.tag         = element_text(size = FONT + 12, face = "bold"),
        plot.margin      = margin(10, 14, 10, 14))

# ── Panel builder ────────────────────────────────────────────────────────────
make_panel <- function(d, title, tag, caption) {
  n_lab <- d %>% count(strat, ejc) %>%
    mutate(label = prettyNum(n, big.mark = ","))

  st <- d %>% group_by(strat) %>%
    wilcox_test(ALLELE.RAT ~ ejc) %>%
    add_significance("p") %>%
    add_xy_position(x = "ejc") %>%
    mutate(y.position = 1.10, ejc = group1)   # ejc column needed by inherited aes

  cat(sprintf("\n[%s] %s\n", tag, title))
  print(n_lab %>% mutate(strat = gsub("\n", " ", strat)) %>%
          select(strat, ejc, n) %>% pivot_wider(names_from = ejc, values_from = n))
  print(st %>% mutate(strat = gsub("\n", " ", strat)) %>% select(strat, n1, n2, p, p.signif))

  ggplot(d, aes(x = ejc, y = ALLELE.RAT, fill = ejc, color = ejc)) +
    geom_jitter(width = 0.18, height = 0, size = 0.9, alpha = 0.25, shape = 16) +
    geom_boxplot(width = 0.38, outlier.shape = NA, alpha = 0.75,
                 color = "grey10", linewidth = 0.8) +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 4.5,
                 fill = "white", color = "grey10", stroke = 1.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50",
               linewidth = 0.6) +
    geom_text(data = n_lab, aes(x = ejc, y = 1.27, label = label),
              inherit.aes = FALSE, size = 6, fontface = "bold", color = "grey10") +
    stat_pvalue_manual(st, label = "p.signif", tip.length = 0.012,
                       bracket.size = 0.7, size = 7, color = "grey20") +
    facet_wrap(~ strat, nrow = 1) +
    scale_fill_manual(values = PAL_EJC) +
    scale_color_manual(values = PAL_EJC) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(-0.02, 1.33),
                       expand = c(0, 0)) +
    labs(y = "NMD efficiency", title = title, tag = tag, caption = caption) +
    base_theme
}

CAP_TX  <- "Transcript-level: median ALLELE.RAT per transcript"
CAP_VAR <- "Variant-level: one point per variant"

pA <- make_panel(to_tx(df,  gene_tol), "Gene tolerance x EJC",   "A", CAP_TX)
pB <- make_panel(to_var(df, af_grp),   "Allele frequency x EJC", "B", CAP_VAR)
pC <- make_panel(to_tx(df,  cds_grp),  "CDS length x EJC",       "C", CAP_TX)
pD <- make_panel(to_tx(df,  exon_grp), "Exon count x EJC",       "D", CAP_TX)

final <- (pA | pB) / (pC | pD)

ggsave("SupplementalFigure3_revised.pdf", final, width = 22, height = 15,
       limitsize = FALSE)
ggsave("SupplementalFigure3_revised.png", final, width = 22, height = 15,
       dpi = 300, limitsize = FALSE)
cat("\nSaved: SupplementalFigure3_revised.pdf / .png\n")
