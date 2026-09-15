#
df_penu <- df_penu %>%
  mutate(
    threeUTR_group = case_when(
      threeUTR.introns == "There is a 3UTR intron" ~ "3UTR-Yes",
      is.na(threeUTR.introns) ~ NA_character_,
      trimws(threeUTR.introns) == "" ~ "3UTR-No",
      TRUE ~ NA_character_
    ),
    threeUTR_group = factor(
      threeUTR_group,
      levels = c("3UTR-No", "3UTR-Yes")
    )
  )

table(df_penu$threeUTR_group, useNA = "ifany")
