# ==============================================================================
# Reviewer #1 Comment #3
# Stratified analysis of PTC-to-EJC distance by 3'UTR intron status
# Penultimate-exon variants
# ==============================================================================
# ==============================================================================
# Reviewer #1 Comment #3
# Stratified analysis of PTC-to-EJC distance by 3'UTR intron status
# Penultimate-exon variants
#
# Purpose:
# Evaluate whether the relationship between PTC-to-EJC distance and NMD
# efficiency differs according to the presence of a 3'UTR intron.
#
# Analysis population:
# - Penultimate-exon variants
# - exon_count >= 3
# - ALLELE.RAT >= 0.35
# - PTC.2.EJC <= 250 nt
#
# Corrected threeUTR.introns coding:
#   "3UTR-No"                 = no 3'UTR intron
#   "There is a 3UTR intron" = 3'UTR intron present
#
# Current filtered counts:
#   No 3'UTR intron      = 520
#   3'UTR intron present = 25
#   Total                 = 545
#
# NOTE:
# LOESS minima are descriptive features of the fitted curve and should not
# be interpreted as formal mechanistic thresholds.
# ==============================================================================


library(tidyverse)
library(ggplot2)
library(mgcv)


# ==============================================================================
# HELPER FUNCTION
# Identify local minima in a LOESS-smoothed curve
# ==============================================================================

find_local_minima <- function(x, y, span = 0.75, n_seq = 400) {

  lo <- loess(
    y ~ x,
    span = span
  )

  x_seq <- seq(
    min(x),
    max(x),
    length.out = n_seq
  )

  y_hat <- predict(
    lo,
    newdata = data.frame(x = x_seq)
  )

  is_min <- c(
    FALSE,
    y_hat[-c(1, n_seq)] < y_hat[-c(n_seq - 1, n_seq)] &
      y_hat[-c(1, n_seq)] < y_hat[-c(1, 2)],
    FALSE
  )

  tibble(
    x_min = x_seq[is_min],
    y_min = y_hat[is_min]
  )
}


# ==============================================================================
# 1. Recode 3'UTR intron status
# ==============================================================================

df_penu <- df_penu %>%
  mutate(

    threeUTR_group = case_when(

      threeUTR.introns == "3UTR-No" ~
        "3UTR-No",

      threeUTR.introns == "There is a 3UTR intron" ~
        "3UTR-Yes",

      TRUE ~ NA_character_
    ),

    threeUTR_group = factor(
      threeUTR_group,
      levels = c(
        "3UTR-No",
        "3UTR-Yes"
      )
    )
  )


cat("\n3'UTR intron status in all penultimate-exon variants:\n")

print(
  table(
    df_penu$threeUTR_group,
    useNA = "always"
  )
)


# ==============================================================================
# 2. Restrict to same PTC-to-EJC range used in Figure 3D
# ==============================================================================

X_MAX_D <- 250


df_D_3UTR <- df_penu %>%
  filter(
    !is.na(PTC.2.EJC),
    !is.na(ALLELE.RAT),
    !is.na(threeUTR_group),
    PTC.2.EJC <= X_MAX_D
  )


df_D_no <- df_D_3UTR %>%
  filter(
    threeUTR_group == "3UTR-No"
  )


df_D_yes <- df_D_3UTR %>%
  filter(
    threeUTR_group == "3UTR-Yes"
  )


cat("\nSample sizes after PTC-to-EJC <= 250 nt filter:\n")

cat(
  "No 3'UTR intron:      ",
  nrow(df_D_no),
  "\n"
)

cat(
  "3'UTR intron present: ",
  nrow(df_D_yes),
  "\n"
)

cat(
  "Total:                 ",
  nrow(df_D_3UTR),
  "\n"
)


# Expected current output:
# No 3'UTR intron:       520
# 3'UTR intron present:   25
# Total:                  545


# ==============================================================================
# 3. Spearman correlations
# ==============================================================================

rho_no <- cor(
  df_D_no$PTC.2.EJC,
  df_D_no$ALLELE.RAT,
  use = "complete.obs",
  method = "spearman"
)


rho_yes <- cor(
  df_D_yes$PTC.2.EJC,
  df_D_yes$ALLELE.RAT,
  use = "complete.obs",
  method = "spearman"
)


cat("\nSpearman correlations:\n")

cat(
  "No 3'UTR intron:      rho =",
  round(rho_no, 3),
  "\n"
)

cat(
  "3'UTR intron present: rho =",
  round(rho_yes, 3),
  "\n"
)


# Current results:
# No 3'UTR intron:       rho =  0.390
# 3'UTR intron present:  rho = -0.165


# ==============================================================================
# 4. LOESS local minima
# ==============================================================================

minima_no <- find_local_minima(
  df_D_no$PTC.2.EJC,
  df_D_no$ALLELE.RAT,
  span = 0.75
)


minima_yes <- find_local_minima(
  df_D_yes$PTC.2.EJC,
  df_D_yes$ALLELE.RAT,
  span = 0.75
)


cat("\nLOESS local minima -- No 3'UTR intron:\n")
print(minima_no)


cat("\nLOESS local minima -- 3'UTR intron present:\n")
print(minima_yes)


# Current results:
#
# No 3'UTR intron:
#   x_min = ~19 nt
#
# 3'UTR intron present:
#   x_min = ~77 nt
#
# Because the 3UTR-Yes subgroup contains only 25 observations, the LOESS
# minimum in that group should be interpreted cautiously and not treated
# as a mechanistic threshold.


# ==============================================================================
# 5. Keep only short-distance minima for figure annotation
# ==============================================================================

minima_no_plot <- minima_no %>%
  filter(
    x_min <= 60
  )


minima_yes_plot <- minima_yes %>%
  filter(
    x_min <= 60
  )


# In the current dataset:
# - No 3'UTR intron: ~19 nt minimum is retained
# - 3'UTR-Yes: ~77 nt minimum is not annotated


# ==============================================================================
# 6. Plot settings
# ==============================================================================

TITLE_COL <- "#2E6DA4"
MIN_COL   <- "#C0392B"
FONT      <- 16


scatter_theme <- theme_classic(
  base_size = FONT
) +
  theme(

    axis.text.x = element_text(
      size = FONT,
      face = "bold",
      color = "grey10"
    ),

    axis.text.y = element_text(
      size = FONT,
      face = "bold",
      color = "grey15"
    ),

    axis.title.x = element_text(
      size = FONT + 3,
      face = "bold",
      margin = margin(t = 12)
    ),

    axis.title.y = element_text(
      size = FONT + 3,
      face = "bold",
      margin = margin(r = 14)
    ),

    panel.grid.major = element_line(
      color = "grey91",
      linewidth = 0.35
    ),

    panel.grid.minor = element_blank(),

    plot.title = element_text(
      size = FONT + 3,
      face = "bold",
      color = TITLE_COL,
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      size = FONT,
      face = "bold.italic",
      color = "grey40",
      hjust = 0.5
    )
  )


y_nmd_scat <- scale_y_continuous(
  breaks = c(
    0,
    0.25,
    0.50,
    0.75,
    1.00
  ),
  limits = c(
    0,
    1.05
  ),
  expand = c(
    0,
    0
  )
)


# ==============================================================================
# 7. Plotting function
# ==============================================================================

make_3UTR_plot <- function(
    dat,
    minima_plot,
    rho,
    subtitle_text
) {

  ggplot(
    dat,
    aes(
      x = PTC.2.EJC,
      y = ALLELE.RAT
    )
  ) +

    geom_point(
      aes(
        color = PTC.2.EJC
      ),
      size = 2.5,
      alpha = 0.45,
      shape = 16
    ) +

    geom_smooth(
      method = "loess",
      span = 0.75,
      color = "black",
      linewidth = 2.0,
      fill = "grey70",
      alpha = 0.25,
      se = TRUE
    ) +

    {
      if (nrow(minima_plot) > 0) list(

        geom_segment(
          data = minima_plot,
          aes(
            x = x_min,
            xend = x_min,
            y = 0,
            yend = y_min
          ),
          color = MIN_COL,
          linetype = "dashed",
          linewidth = 1.2,
          inherit.aes = FALSE
        ),

        geom_point(
          data = minima_plot,
          aes(
            x = x_min,
            y = y_min
          ),
          shape = 25,
          size = 4,
          fill = MIN_COL,
          color = "white",
          stroke = 1.2,
          inherit.aes = FALSE
        ),

        geom_text(
          data = minima_plot,
          aes(
            x = x_min + 6,
            y = 0.32,
            label = sprintf(
              "%.0f nt",
              x_min
            )
          ),
          hjust = 0,
          size = 5,
          fontface = "bold.italic",
          color = MIN_COL,
          inherit.aes = FALSE
        )
      )
    } +

    geom_hline(
      yintercept = 0.5,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.8
    ) +

    annotate(
      "text",
      x = 155,
      y = 0.08,
      label = sprintf(
        "\u03c1 = %.3f",
        rho
      ),
      size = 5,
      fontface = "bold.italic",
      color = "grey20"
    ) +

    annotate(
      "text",
      x = 3,
      y = 0.98,
      label = paste0(
        "n = ",
        nrow(dat)
      ),
      hjust = 0,
      size = 5,
      fontface = "bold",
      color = "grey20"
    ) +

    scale_color_gradient2(
      low = "#2980B9",
      mid = "#F0E442",
      high = "#C0392B",
      midpoint = 125,
      name = "PTC-to-EJC\n(nt)",
      guide = guide_colorbar(
        barwidth = 1.3,
        barheight = 9,
        title.hjust = 0.5
      )
    ) +

    scale_x_continuous(
      breaks = seq(
        0,
        250,
        50
      ),
      limits = c(
        -2,
        258
      ),
      expand = c(
        0.01,
        0
      )
    ) +

    y_nmd_scat +

    labs(
      x = "PTC-to-EJC distance (nt)",
      y = "NMD efficiency",
      title = "PTC-to-EJC Distance vs NMD Efficiency",
      subtitle = subtitle_text
    ) +

    scatter_theme +

    theme(
      legend.position = "right",

      legend.title = element_text(
        size = FONT - 1,
        face = "bold"
      ),

      legend.text = element_text(
        size = FONT - 2
      )
    )
}


# ==============================================================================
# 8. Generate stratified figures
# ==============================================================================

p_no3UTR <- make_3UTR_plot(
  dat = df_D_no,
  minima_plot = minima_no_plot,
  rho = rho_no,
  subtitle_text = "No 3\u2032UTR intron"
)


p_yes3UTR <- make_3UTR_plot(
  dat = df_D_yes,
  minima_plot = minima_yes_plot,
  rho = rho_yes,
  subtitle_text = "3\u2032UTR intron present"
)


print(p_no3UTR)
print(p_yes3UTR)


# ==============================================================================
# 9. Save figures
# ==============================================================================

ggsave(
  "PTC_EJC_penultimate_No_3UTR_intron.pdf",
  p_no3UTR,
  width = 10,
  height = 8,
  units = "in"
)


ggsave(
  "PTC_EJC_penultimate_3UTR_intron_present.pdf",
  p_yes3UTR,
  width = 10,
  height = 8,
  units = "in"
)


# ==============================================================================
# 10. GAM analysis
# Test whether 3'UTR-intron status modifies the distance-response relationship
# ==============================================================================

df_gam <- df_D_3UTR %>%
  mutate(

    # Treatment-coded factor for interpretable main group effect
    threeUTR_fac = factor(
      threeUTR_group,
      levels = c(
        "3UTR-No",
        "3UTR-Yes"
      )
    ),

    # Ordered factor used for the difference smooth
    threeUTR_ord = ordered(
      threeUTR_group,
      levels = c(
        "3UTR-No",
        "3UTR-Yes"
      )
    )
  )


# ------------------------------------------------------------------------------
# Shared-curve model
# Assumes same PTC-to-EJC distance-response relationship in both groups
# ------------------------------------------------------------------------------

m_same <- gam(
  ALLELE.RAT ~
    threeUTR_fac +
    s(
      PTC.2.EJC,
      k = 5
    ),
  data = df_gam,
  method = "REML"
)


# ------------------------------------------------------------------------------
# Difference-smooth model
# Allows 3UTR-Yes variants to have a different distance-response relationship
# ------------------------------------------------------------------------------

m_diff2 <- gam(
  ALLELE.RAT ~
    threeUTR_fac +
    s(
      PTC.2.EJC,
      k = 5
    ) +
    s(
      PTC.2.EJC,
      by = threeUTR_ord,
      k = 5
    ),
  data = df_gam,
  method = "REML"
)


cat("\n============================================================\n")
cat("GAM: SHARED-CURVE MODEL\n")
cat("============================================================\n")

print(
  summary(m_same)
)


cat("\n============================================================\n")
cat("GAM: DIFFERENCE-SMOOTH MODEL\n")
cat("============================================================\n")

print(
  summary(m_diff2)
)


# ==============================================================================
# 11. ML sensitivity analysis for model comparison
# ==============================================================================

m_same_ML <- gam(
  ALLELE.RAT ~
    threeUTR_fac +
    s(
      PTC.2.EJC,
      k = 5
    ),
  data = df_gam,
  method = "ML"
)


m_diff_ML <- gam(
  ALLELE.RAT ~
    threeUTR_fac +
    s(
      PTC.2.EJC,
      k = 5
    ) +
    s(
      PTC.2.EJC,
      by = threeUTR_ord,
      k = 5
    ),
  data = df_gam,
  method = "ML"
)


cat("\n============================================================\n")
cat("ML MODEL COMPARISON\n")
cat("============================================================\n")


print(
  anova(
    m_same_ML,
    m_diff_ML,
    test = "F"
  )
)


cat("\nAIC comparison:\n")

print(
  AIC(
    m_same_ML,
    m_diff_ML
  )
)


# ==============================================================================
# 12. GAM diagnostics
# ==============================================================================

cat("\n============================================================\n")
cat("GAM DIAGNOSTICS\n")
cat("============================================================\n")

gam.check(
  m_diff2
)


# ==============================================================================
# CURRENT RESULTS
# ==============================================================================
#
# Sample sizes:
#   No 3'UTR intron      = 520
#   3'UTR intron present = 25
#   Total                 = 545
#
# Spearman correlation:
#   No 3'UTR intron       rho =  0.390
#   3'UTR intron present  rho = -0.165
#
# LOESS local minimum:
#   No 3'UTR intron       ~19 nt
#   3'UTR intron present  ~77 nt
#
# NOTE:
# The ~77 nt minimum in the 3UTR-Yes subgroup should not be interpreted as
# a biological threshold because only 25 observations are available and the
# LOESS confidence interval is wide.
#
#
# GAM difference-smooth model:
#
# Overall PTC-to-EJC distance smooth:
#   edf = 3.544
#   F   = 26.401
#   p   < 2e-16
#
# 3UTR-Yes difference smooth:
#   edf = 1.000
#   F   = 5.272
#   p   = 0.022
#
# Main 3UTR group term:
#   Estimate = 0.0825
#   p = 0.0147
#
# Model:
#   adjusted R-squared = 0.163
#   deviance explained = 17.2%
#
#
# ML nested-model comparison:
#   F = 5.241
#   p = 0.0224
#
# AIC:
#   shared curve        = -412.87
#   group-specific curve = -416.15
#
# Lower AIC favors the model allowing the PTC-to-EJC relationship
# to differ by 3'UTR-intron status.
#
#
# Interpretation:
#
# The approximately 18-19 nt short-distance feature observed in Figure 3D
# persists among the 520 variants without a 3'UTR intron, indicating that
# the presence of a 3'UTR intron is not required to produce this feature.
#
# However, the GAM provides evidence that the overall PTC-to-EJC
# distance-response relationship differs according to 3'UTR-intron status
# (p = 0.022).
#
# Because only 25 variants contain a 3'UTR intron, the subgroup-specific
# result should be interpreted cautiously and requires confirmation in
# larger datasets.
#
# ==============================================================================

