############################################################
# Reviewer #1, Comment #3
#
# Figure 3D: PTC-to-EJC distance stratified by 3'UTR introns
#
# Analyses:
# 1. Stratified LOESS curves and Spearman correlations
# 2. LOESS local minima
# 3. GAM interaction and model comparison
#
# Dataset: TOPMed stop-gain variants
############################################################

library(dplyr)
library(ggplot2)
library(mgcv)
library(patchwork)

# ----------------------------------------------------------
# 1. Load and filter data
# ----------------------------------------------------------

df <- read.csv(
  "TOPMed_stopgain_September25_corrected_readyformodel.csv"
)

df_D <- df %>%
  filter(
    penultimate.full == "penultimate",
    exon_count >= 3,
    ALLELE.RAT >= 0.35,
    !is.na(PTC.2.EJC),
    between(PTC.2.EJC, 0, 250)
  ) %>%
  mutate(
    threeUTR_group = case_when(
      threeUTR.introns == "3UTR-No" ~ "No 3'UTR intron",
      threeUTR.introns == "There is a 3UTR intron" ~
        "3'UTR intron present",
      TRUE ~ NA_character_
    ),
    threeUTR_group = factor(
      threeUTR_group,
      levels = c(
        "3'UTR intron present",
        "No 3'UTR intron"
      )
    )
  ) %>%
  filter(!is.na(threeUTR_group))

print(table(df_D$threeUTR_group))

# Historical expected counts:
# 3'UTR intron present: 25
# No 3'UTR intron:      520


# ----------------------------------------------------------
# 2. Spearman correlations
# ----------------------------------------------------------

stats <- df_D %>%
  group_by(threeUTR_group) %>%
  summarise(
    n = n(),
    rho = cor(
      PTC.2.EJC,
      ALLELE.RAT,
      method = "spearman"
    ),
    .groups = "drop"
  )

print(stats)


# ----------------------------------------------------------
# 3. Find LOESS local minima
# ----------------------------------------------------------

find_minima <- function(dat, span = 0.75) {

  fit <- loess(
    ALLELE.RAT ~ PTC.2.EJC,
    data = dat,
    span = span
  )

  x <- seq(
    min(dat$PTC.2.EJC),
    max(dat$PTC.2.EJC),
    length.out = 400
  )

  y <- predict(
    fit,
    newdata = data.frame(PTC.2.EJC = x)
  )

  idx <- which(
    diff(sign(diff(y))) == 2
  ) + 1

  data.frame(
    PTC.2.EJC = x[idx],
    ALLELE.RAT = y[idx]
  )
}

minima <- bind_rows(
  lapply(
    split(df_D, df_D$threeUTR_group),
    function(dat) {
      result <- find_minima(dat)
      result$threeUTR_group <- as.character(
        dat$threeUTR_group[1]
      )
      result
    }
  )
) %>%
  mutate(
    threeUTR_group = factor(
      threeUTR_group,
      levels = levels(df_D$threeUTR_group)
    )
  )

print(minima)

# Mark only the short-distance minimum, when present.
short_minima <- minima %>%
  filter(between(PTC.2.EJC, 10, 60))


# ----------------------------------------------------------
# 4. Stratified LOESS figure
# ----------------------------------------------------------

p <- ggplot(
  df_D,
  aes(PTC.2.EJC, ALLELE.RAT)
) +
  geom_point(
    aes(color = PTC.2.EJC),
    size = 2,
    alpha = 0.5
  ) +
  geom_smooth(
    method = "loess",
    span = 0.75,
    color = "black",
    fill = "grey75",
    se = TRUE
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_vline(
    data = short_minima,
    aes(xintercept = PTC.2.EJC),
    color = "red3",
    linetype = "dashed"
  ) +
  geom_point(
    data = short_minima,
    aes(PTC.2.EJC, ALLELE.RAT),
    inherit.aes = FALSE,
    color = "red3",
    size = 3
  ) +
  facet_wrap(
    ~threeUTR_group,
    nrow = 1
  ) +
  scale_color_gradient(
    low = "#2980B9",
    high = "#C0392B",
    name = "PTC-to-EJC (nt)"
  ) +
  scale_x_continuous(
    breaks = seq(0, 250, 50)
  ) +
  coord_cartesian(
    xlim = c(0, 250),
    ylim = c(0, 1.05)
  ) +
  labs(
    x = "PTC-to-EJC distance (nt)",
    y = "NMD efficiency"
  ) +
  theme_classic(base_size = 14)

# Add sample sizes and Spearman correlations
p <- p +
  geom_text(
    data = stats,
    aes(
      x = 125,
      y = 0.08,
      label = sprintf("n = %d, rho = %.3f", n, rho)
    ),
    inherit.aes = FALSE,
    size = 4
  )

print(p)

ggsave(
  "Figure3D_3UTR_stratification.pdf",
  p,
  width = 11,
  height = 5
)


# ----------------------------------------------------------
# 5. GAM: Does 3'UTR intron status modify the relationship?
# ----------------------------------------------------------

df_gam <- df_D %>%
  mutate(
    threeUTR_fac = factor(
      threeUTR_group,
      levels = c(
        "No 3'UTR intron",
        "3'UTR intron present"
      )
    ),
    threeUTR_ord = ordered(
      threeUTR_group,
      levels = c(
        "No 3'UTR intron",
        "3'UTR intron present"
      )
    )
  )

# Shared distance-response curve
m_same <- gam(
  ALLELE.RAT ~
    threeUTR_fac +
    s(PTC.2.EJC, k = 5),
  data = df_gam,
  method = "ML"
)

# Additional smooth for 3'UTR-intron group
m_diff <- gam(
  ALLELE.RAT ~
    threeUTR_fac +
    s(PTC.2.EJC, k = 5) +
    s(
      PTC.2.EJC,
      by = threeUTR_ord,
      k = 5
    ),
  data = df_gam,
  method = "ML"
)

cat("\nGAM with group-specific smooth:\n")
print(summary(m_diff))

cat("\nGAM model comparison:\n")
print(anova(m_same, m_diff, test = "F"))

cat("\nAIC comparison:\n")
print(AIC(m_same, m_diff))
