#PTC-2-start
breaks_start <- c(0,100,200,300,700, max(df$PTC.2.start, na.rm = TRUE))
# Define matching labels
labels_start <- c("0-100","100-200","200-300","300-700",">700")
# Apply binning
df$PTC.2.start.binning <- cut(
    df$PTC.2.start,
    breaks = breaks_start,
    labels = labels_start,
    include.lowest = TRUE,
    right = FALSE   
)
#PTC-2-end
breaks_end <- c(0,100,200,300,700, max(df$PTC.2.start, na.rm = TRUE))
# Define matching labels
labels_end <- c("0-100","100-200","200-300","300-700",">700")

# Apply binning
df$PTC.2.end.binning <- cut(
    df$PTC.2.end,
    breaks = breaks_end,
    labels = labels_end,
    include.lowest = TRUE,
    right = FALSE   
)
#PTC-2-EJC
breaks_ejc <- c(0,100,200,300,700, max(df$PTC.2.start, na.rm = TRUE))
# Define matching labels
labels_ejc <- c("0-100","100-200","200-300","300-700",">700")
# Apply binning
df$PTC.2.EJC.binning <- cut(
    df$PTC.2.EJC,
    breaks = breaks_ejc,
    labels = labels_ejc,
    include.lowest = TRUE,
    right = FALSE   
)
#PTC bearing exon
breaks_exon <- c(0,100,200,300,700,max(df$length.mutated.exon, na.rm = TRUE))

# Define matching labels
labels_exon <- c("0-100","100-200","200-300", "300-700",">700")
# Apply binning
df$length.mutated.exon.binning <- cut(
    df$length.mutated.exon,
    breaks = breaks_exon,
    labels = labels_exon,
    include.lowest = TRUE,
    right = FALSE  
)

#PTC (amino acid)
df$minus1_is_G <- rep('Not Glycine',nrow(df))
ind.1 <- which(df$PTC.minus1.aa=='G')
df$minus1_is_G[ind.1] <- rep('Glycine',length(ind.1))


df$minus2_is_G <- rep('Not Glycine',nrow(df))
ind.1 <- which(df$PTC.minus2.aa=='G')
df$minus2_is_G[ind.1] <- rep('Glycine',length(ind.1))

df$plus1_is_G <- rep('Not Glycine',nrow(df))
ind.1 <- which(df$PTC.plus1.aa=='G')
df$plus1_is_G[ind.1] <- rep('Glycine',length(ind.1))

df$plus2_is_G <- rep('Not Glycine',nrow(df))
ind.1 <- which(df$PTC.plus2.aa=='G')
df$plus2_is_G[ind.1] <- rep('Glycine',length(ind.1))


#Dataset uploading

new_df2 <- jake1   # start once

for (col in phastCons_cols) {
  new_col <- paste0(col, "_cat")
  
  new_df2 <- categorise_score(
    new_df2,   # ← NOT jake1 anymore
    col = col,
    breaks = c(0.1, 0.5),
    labels = c("Low (<0.1)", "Medium (0.1=< & <0.5)", "High (>=0.5)"),
    new_col_name = new_col
  )
}


for (col in phyloP_cols) {
  new_col <- paste0(col, "_cat")
  
  new_df2 <- categorise_score(
    new_df2,
    col = col,
    breaks = c(0.5, 2),
    labels = c("Low/Neutral (<0.5)", "Medium (0.5=< & <2)", "High (>=2)"),
    new_col_name = new_col
  )
}

table(new_df2$phastcons_ejc_100bp_median_cat)






