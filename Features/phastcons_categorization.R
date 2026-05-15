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




