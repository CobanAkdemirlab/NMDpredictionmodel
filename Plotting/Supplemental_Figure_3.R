# ══════════════════════════════════════════════════════════════════════════════
# NMD Supplemental Figure 3 — 4-panel transcript features (v3)
#
# v3: unit of analysis restored to match the original (June) Figure S2:
#   A, C, D = TRANSCRIPT-level: one point per transcript, median ALLELE.RAT of
#             its variants (original n: 348 / 2,739; 1,301 / 1,028 / 883;
#             694 / 1,024 / 1,021 / 473 -> 3,212 transcripts)
#   B       = VARIANT-level: allele frequency is a property of the variant, not
#             the transcript (original n: 4,163 / 1,584)
#   v2 had plotted every panel at variant level. Each panel now carries a
#   caption stating its level, as in Supplemental Figure 3.
#
# v2 (corrected dataset):
#   - Input switched from TOPMed_stopgain_withfullpenultimateMarch30.csv to the
#     corrected TOPMed_stopgain_September25_corrected_readyformodel.csv
#     (same file as revised Figures 3, 4 and 6)
#   - Filter changed from ALLELE.RAT > 0.35 to ALLELE.RAT >= 0.35, matching the
#     main figures (5,749 variants)
#   - Checks that every required column exists before plotting; if pLI.cat is
#     missing but a numeric pLI column is present, the tolerance bins are
#     rebuilt from it with the same cut-offs (>= 0.65 / < 0.35)
#   - Prints the n per group for each panel (for the figure legend)
#   - Output renamed SupplementalFigure2_revised.pdf / .png
#   - Layout: n labels raised above the bracket tiers (no collisions), no
#     padded spaces in n labels, exon-count labels on two lines
#
# Panels (single row):
#   A = Gene tolerance  (Highly intolerant vs Highly tolerant)        [Mann-Whitney]
#   B = Allele frequency (Ultra-rare vs Rare/Common)                   [Mann-Whitney]
#   C = CDS length      (Short / Medium / Long)                        [Kruskal-Wallis]
#   D = Exon count      (2-5 / 6-10 / 11-20 / >20)                     [Kruskal-Wallis]
#   All panels are variant-level: one point per variant, no per-transcript
#   aggregation.
#
# Required file:
#   TOPMed_stopgain_September25_corrected_readyformodel.csv
#
# Filter applied:
#   ALLELE.RAT >= 0.35
#
# Expected columns (as in the March 30 file; checked at load time):
#   pLI.cat     -> already-binned gene tolerance
#                  ("highly intolerant (pLI>=0.65)" / "highly tolerant (pLI<0.35)")
#   Freq.cat    -> already-binned allele frequency
#                  ("Ultra-rare variants" / "Rare/Common variants")
#   cds_length  -> numeric, in bp; cuts at 1289 / 2181
#   exon_count  -> numeric, total exons; cuts at 5 / 10 / 20
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(ggplot2)
library(patchwork)
library(rstatix)
library(ggpubr)

# ══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ══════════════════════════════════════════════════════════════════════════════
TITLE_COL  <- "#2E6DA4"
FONT       <- 22
PANEL_FILL <- "#7BB6A0"
BOX_FILL   <- "white"

CDS_CUTS   <- c(1289, 2181)
EXON_CUTS  <- c(5, 10, 20)

# ══════════════════════════════════════════════════════════════════════════════
# LOAD & FILTER
# ══════════════════════════════════════════════════════════════════════════════
INPUT_CSV <- "TOPMed_stopgain_September25_corrected_readyformodel.csv"
AR_MIN    <- 0.35

df <- read.csv(INPUT_CSV, stringsAsFactors = FALSE) %>%
  filter(!is.na(ALLELE.RAT),
         ALLELE.RAT >= AR_MIN)          # >= (was >), same as Figures 3, 4, 6

cat(sprintf("\n[Input: %s]\n", INPUT_CSV))
cat(sprintf("[Filter applied: ALLELE.RAT >= %.2f]\n", AR_MIN))
cat(sprintf("Filtered variants: %d (main figures use 5,749)\n", nrow(df)))

# ── Column check ─────────────────────────────────────────────────────────────
# Gene tolerance: use pLI.cat if present; otherwise rebuild it from numeric pLI
if (!"pLI.cat" %in% names(df)) {
  pli_col <- intersect(c("pLI", "pli", "gnomad_pLI", "pLI_score"), names(df))
  if (length(pli_col) == 0)
    stop("Neither 'pLI.cat' nor a numeric pLI column was found. Columns containing ",
         "'pli': ", paste(grep("pli", names(df), ignore.case = TRUE, value = TRUE),
                          collapse = ", "))
  cat(sprintf("pLI.cat not found; rebuilding tolerance bins from '%s'\n", pli_col[1]))
  df$pLI.cat <- case_when(
    df[[pli_col[1]]] >= 0.65 ~ "highly intolerant (pLI>=0.65)",
    df[[pli_col[1]]] <  0.35 ~ "highly tolerant (pLI<0.35)",
    TRUE ~ NA_character_)
}

required <- c("pLI.cat", "Freq.cat", "cds_length", "exon_count")
missing  <- setdiff(required, names(df))
if (length(missing) > 0)
  stop("Column(s) not found in ", INPUT_CSV, ": ", paste(missing, collapse = ", "),
       "\nPossible matches: ",
       paste(grep("freq|af|pli|cds|exon", names(df), ignore.case = TRUE, value = TRUE),
             collapse = ", "))

# Category labels must match exactly, or the group silently becomes NA
cat("\npLI.cat values:\n");  print(table(df$pLI.cat,  useNA = "ifany"))
cat("Freq.cat values:\n");   print(table(df$Freq.cat, useNA = "ifany"))

# ══════════════════════════════════════════════════════════════════════════════
# DERIVE GROUPING VARIABLES (using verified column names)
# ══════════════════════════════════════════════════════════════════════════════
df <- df %>%
  mutate(
    # A — Gene tolerance: pLI.cat is already a binned string
    gene_tol = case_when(
      pLI.cat == "highly intolerant (pLI>=0.65)" ~ "Highly\nintolerant",
      pLI.cat == "highly tolerant (pLI<0.35)"    ~ "Highly\ntolerant",
      TRUE                                        ~ NA_character_
    ),
    gene_tol = factor(gene_tol,
                      levels = c("Highly\nintolerant","Highly\ntolerant")),

    # B — Allele frequency: Freq.cat is already a binned string
    af_grp = case_when(
      Freq.cat == "Ultra-rare variants"  ~ "Ultra-rare",
      Freq.cat == "Rare/Common variants" ~ "Rare/Common",
      TRUE                                ~ NA_character_
    ),
    af_grp = factor(af_grp, levels = c("Ultra-rare","Rare/Common")),

    # C — CDS length, cuts at 1289 / 2181 bp
    cds_grp = case_when(
      cds_length <  CDS_CUTS[1] ~ "Short\n<1.3 kb",
      cds_length <= CDS_CUTS[2] ~ "Medium\n1.3-2.2 kb",
      TRUE                      ~ "Long\n>2.2 kb"
    ),
    cds_grp = factor(cds_grp,
                     levels = c("Short\n<1.3 kb",
                                "Medium\n1.3-2.2 kb",
                                "Long\n>2.2 kb")),

    # D — Exon count, cuts at 5 / 10 / 20
    exon_grp = case_when(
      exon_count <= EXON_CUTS[1] ~ "2-5\nexons",
      exon_count <= EXON_CUTS[2] ~ "6-10\nexons",
      exon_count <= EXON_CUTS[3] ~ "11-20\nexons",
      TRUE                       ~ ">20\nexons"
    ),
    exon_grp = factor(exon_grp,
                      levels = c("2-5\nexons","6-10\nexons",
                                 "11-20\nexons",">20\nexons"))
  )

stopifnot("No variants matched the pLI.cat labels"  = any(!is.na(df$gene_tol)),
          "No variants matched the Freq.cat labels" = any(!is.na(df$af_grp)))

# ── v3: transcript-level table (one row per transcript) for panels A, C, D ──
# gene tolerance, CDS length and exon count are transcript properties
tx <- df %>%
  group_by(TxName) %>%
  summarise(ALLELE.RAT = median(ALLELE.RAT),
            gene_tol   = first(gene_tol),
            cds_grp    = first(cds_grp),
            exon_grp   = first(exon_grp),
            n_var      = n(), .groups = "drop")
cat(sprintf("\nTranscripts: %d  (original Figure S2: 3,212)\n", nrow(tx)))

# n per group (for the Supplemental Figure 2 legend)
cat("\n── n per group (for the legend) ──\n")
for (v in c("gene_tol", "af_grp", "cds_grp", "exon_grp")) {
  src <- if (v == "af_grp") df else tx
  cat(v, if (v == "af_grp") "(variant-level):\n" else "(transcript-level):\n")
  print(src %>% count(.data[[v]]) %>%
          mutate(group = gsub("\n", " ", as.character(.data[[v]]))) %>%
          select(group, n), row.names = FALSE)
}

# ══════════════════════════════════════════════════════════════════════════════
# THEME
# ══════════════════════════════════════════════════════════════════════════════
base_theme <- theme_classic(base_size = FONT) +
  theme(
    legend.position = "none",
    plot.title      = element_text(size = FONT + 2, face = "bold",
                                   color = TITLE_COL, hjust = 0.5,
                                   lineheight = 1.1, margin = margin(b = 8)),
    axis.title.y    = element_text(size = FONT, face = "bold"),
    axis.title.x    = element_blank(),
    axis.text.x     = element_text(size = FONT - 2, face = "bold",
                                   color = "grey10", lineheight = 1.0),
    axis.text.y     = element_text(size = FONT - 2, color = "grey10"),
    axis.line       = element_line(linewidth = 0.9, color = "grey20"),
    axis.ticks      = element_line(linewidth = 0.9, color = "grey20"),
    plot.tag        = element_text(size = FONT + 14, face = "bold"),
    plot.margin     = margin(8, 14, 8, 14)
  )

ref05 <- geom_hline(yintercept = 0.5, linetype = "dashed",
                    color = "grey55", linewidth = 0.7)

# ══════════════════════════════════════════════════════════════════════════════
# HELPER — build a violin+box panel
# ══════════════════════════════════════════════════════════════════════════════
make_panel <- function(data, x, panel_title, panel_tag,
                       comparisons = NULL, kruskal = FALSE, caption = NULL) {
  data <- data %>% drop_na({{ x }})
  n_lbls <- data %>% count({{ x }}) %>%
    rename(label_grp = 1, n = 2) %>%
    mutate(label = paste0("n = ", prettyNum(n, big.mark = ",")))   # v2: no padding spaces

  if (kruskal) {
    p <- kruskal_test(data, as.formula(paste("ALLELE.RAT ~",
                                              rlang::as_name(rlang::enquo(x)))))$p
    test_lbl <- sprintf("Kruskal-Wallis: %s",
                        if (p < 0.0001) "****" else if (p < 0.001) "***"
                        else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns")
  } else {
    p <- wilcox_test(data, as.formula(paste("ALLELE.RAT ~",
                                             rlang::as_name(rlang::enquo(x)))))$p
    test_lbl <- sprintf("Mann-Whitney: %s",
                        if (p < 0.0001) "****" else if (p < 0.001) "***"
                        else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns")
  }

  stat_df <- NULL
  if (!is.null(comparisons)) {
    stat_df <- data %>%
      wilcox_test(as.formula(paste("ALLELE.RAT ~",
                                    rlang::as_name(rlang::enquo(x)))),
                  comparisons = comparisons,
                  p.adjust.method = "BH") %>%
      add_significance("p") %>%
      add_xy_position(x = rlang::as_name(rlang::enquo(x)))
    if (nrow(stat_df) > 1) {
      # v2: two tiers (1.05 / 1.13) kept below the n labels at 1.26
      stat_df$y.position <- 1.05 + 0.08 * (seq_len(nrow(stat_df)) - 1) %% 2
    } else {
      stat_df$y.position <- 1.08
    }
  }

  p_full <- ggplot(data, aes(x = {{ x }}, y = ALLELE.RAT)) +
    geom_violin(fill = PANEL_FILL, color = "grey10",
                alpha = 0.78, linewidth = 1.2,
                trim = TRUE, width = 0.85) +
    geom_boxplot(width = 0.22, outlier.shape = NA,
                 fill = BOX_FILL, color = "grey10",
                 linewidth = 1.2, alpha = 0.95) +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 4.5,
                 fill = "white", color = "grey10", stroke = 1.4) +
    ref05 +
    geom_text(data = n_lbls, aes(x = label_grp, y = 1.26, label = label),
              inherit.aes = FALSE, size = 7, fontface = "bold",
              color = "grey10")

  if (!is.null(stat_df)) {
    p_full <- p_full +
      stat_pvalue_manual(stat_df, label = "p.signif",
                         tip.length = 0.012, bracket.size = 1.0,
                         size = 10, color = "grey20")
  }

  p_full +
    scale_y_continuous(breaks = seq(0, 1, 0.25),
                       limits = c(-0.04, 1.32), expand = c(0, 0)) +
    labs(y = "NMD efficiency",
         title = paste0(panel_title, "\n(", test_lbl, ")"),
         tag = panel_tag, caption = caption) +
    base_theme +
    theme(plot.caption = element_text(size = FONT - 7, color = "grey45",
                                      hjust = 0.5))
}

# ══════════════════════════════════════════════════════════════════════════════
# BUILD PANELS
# ══════════════════════════════════════════════════════════════════════════════
pA <- make_panel(tx, gene_tol, panel_title = "Gene tolerance", panel_tag = "A",
                 comparisons = list(c("Highly\nintolerant","Highly\ntolerant")),
                 kruskal = FALSE, caption = "Transcript-level: median ALLELE.RAT per transcript")

pB <- make_panel(df, af_grp, panel_title = "Allele frequency", panel_tag = "B",
                 comparisons = list(c("Ultra-rare","Rare/Common")),
                 kruskal = FALSE, caption = "Variant-level: one point per variant")

pC <- make_panel(tx, cds_grp, panel_title = "CDS length", panel_tag = "C",
                 comparisons = list(c("Short\n<1.3 kb","Medium\n1.3-2.2 kb"),
                                     c("Medium\n1.3-2.2 kb","Long\n>2.2 kb")),
                 kruskal = TRUE, caption = "Transcript-level: median ALLELE.RAT per transcript")

pD <- make_panel(tx, exon_grp, panel_title = "Exon count", panel_tag = "D",
                 comparisons = list(c("2-5\nexons","6-10\nexons"),
                                     c("6-10\nexons","11-20\nexons"),
                                     c("11-20\nexons",">20\nexons")),
                 kruskal = TRUE, caption = "Transcript-level: median ALLELE.RAT per transcript")

# ══════════════════════════════════════════════════════════════════════════════
# COMBINE & SAVE
# ══════════════════════════════════════════════════════════════════════════════
final <- (pA | pB | pC | pD) &
  theme(plot.margin = margin(10, 14, 10, 14))

ggsave("SupplementalFigure2_revised.pdf", final,
       width = 32, height = 9, dpi = 300, limitsize = FALSE)
ggsave("SupplementalFigure2_revised.png", final,
       width = 32, height = 9, dpi = 300, limitsize = FALSE)

print(final)
cat("\nSaved: SupplementalFigure2_revised.pdf / .png\n")
