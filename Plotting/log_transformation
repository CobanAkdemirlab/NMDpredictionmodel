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
