# ==============================================================================
# Exploratory feature transformations
#
# Purpose:
#   Generate log2-transformed and categorized versions of selected transcript
#   length features for descriptive analyses, single-feature association tests,
#   and plotting.
#
# These variables were not used as primary transcript annotations and should
# not be confused with the continuous features used in the feature-generation
# workflow.
#
# Input variables:
#   - cds_length
#   - threeUTR_length
#   - fiveutr_length
#   - log2newUTR
#   - pLI
#
# Derived exploratory variables:
#   - log2_CDS
#   - cds_length.cut
#   - log2_3utr
#   - threeUTR_length.cut
#   - log2_5utr
#   - fiveUTR_length.cut
#   - log2newUTR.cut
#   - pLI.cat
# ==============================================================================
df$log2_CDS <- log2(df$cds_length)
max_log2 <- max(df$log2_CDS, na.rm = TRUE)

df$cds_length.cut <- cut(
  df$log2_CDS,
  breaks = c(0, 10, max_log2),
  labels = c("0-10", ">10"),
  include.lowest = TRUE,
  right = TRUE
)

# Check distribution
table(df$cds_length.cut, useNA = "ifany")


df$log2_3utr <- log2(df$threeUTR_length)
max_log2__3utr <- max(df$log2_3utr, na.rm = TRUE)

df$threeUTR_length.cut<- cut(
  df$log2_3utr,
  breaks = c(0, 10, max_log2__3utr),
  labels = c("0-10", ">10"),
  include.lowest = TRUE,
  right = TRUE
)

# Check distribution
table(df$threeUTR_length.cut, useNA = "ifany")

df$log2_5utr <- log2(df$fiveutr_length)
max_log2__5utr <- max(df$log2_5utr, na.rm = TRUE)

df$fiveUTR_length.cut<- cut(
  df$log2_5utr,
  breaks = c(0, 7, max_log2__5utr),
  labels = c("0-7", ">7"),
  include.lowest = TRUE,
  right = TRUE
)

# Check distribution
table(df$fiveUTR_length.cut, useNA = "ifany")

#New 3'UTR 
#df$newUTR_length <- df$threeUTR_length+df$PTC.2.end
#df$log2newUTR <- log2(df$newUTR_length)

max_log2newUTR <- max(df$log2newUTR, na.rm = TRUE)
df$log2newUTR.cut <- cut(
  df$log2newUTR,
  breaks = c(0, 10, max_log2newUTR),
  labels = c("0-10", ">10"),
  include.lowest = TRUE,
  right = TRUE
)

# Check distribution
table(df$log2newUTR.cut, useNA = "ifany")

# ------------------------------------------------------------------------------
# pLI categories
# Used for exploratory analyses and visualization only.
# The continuous pLI score is retained as the gene-level model feature.
# ------------------------------------------------------------------------------

df$pLI.cat <- NA_character_

df$pLI.cat[df$pLI < 0.35] <-
  "highly tolerant (pLI < 0.35)"

df$pLI.cat[df$pLI >= 0.35 & df$pLI < 0.65] <-
  "medium tolerant (0.35 <= pLI < 0.65)"

df$pLI.cat[df$pLI >= 0.65] <-
  "highly intolerant (pLI >= 0.65)"


