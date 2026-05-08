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


#dataset uploading


# A general categorise function
categorise_score <- function(df, col, breaks, labels, new_col_name) {
  jake1 %>%
    mutate(
      !!new_col_name := case_when(
        (!!sym(col) < breaks[1])                               ~ labels[1],
        (!!sym(col) >= breaks[1] & !!sym(col) < breaks[2])     ~ labels[2],
        (!!sym(col) >= breaks[2])                              ~ labels[3],
        TRUE                                                    ~ NA_character_
      )
    )
}


phastCons_cols <- c("phastcons_ejc_100bp_median",
                    "phastcons_new3utr_first200_median",
                    "phastcons_new3utr_whole_median",
                    "phastcons_old3utr_first200_median",
                    "phastcons_old3utr_whole_median",
                    "phastcons_ptc_100bp_median",
                    "phastcons_ptc_to_ejc_median",
                    "phastcons_tx_whole_median",
                    "phastcons_utr5_first200_median",
                    "phastcons_utr5_whole_median")

for (col in phastCons_cols) {
  new_col <- paste0(col, "_cat")
  new_df2 <- categorise_score(jake1, col = col,
                              breaks = c(0.1, 0.5),
                              labels = c("Low", "Medium", "High"),
                              new_col_name = new_col)

}


phyloP_cols <- c("phylop_new3utr_first200_median",
                 "phylop_new3utr_whole_median",
                 "phylop_old3utr_first200_median",
                 "phylop_old3utr_whole_median",
                 "phylop_ptc_100bp_median",
                 "phylop_ptc_to_ejc_median",
                 "phylop_tx_whole_median",
                 "phylop_utr5_first200_median",
                 "phylop_ejc_100bp_median",
                 "phylop_utr5_whole_median")

for (col in phyloP_cols) {
  new_col <- paste0(col, "_cat")
  new_df2 <- categorise_score(jake1, col = col,
                              breaks = c(0.5, 2),
                              labels = c("Low/Neutral", "Medium", "High"),
                              new_col_name = new_col)
}




