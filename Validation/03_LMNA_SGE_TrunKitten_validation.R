############################################################
# LMNA SGE experimental validation of TrunKitten
############################################################

library(dplyr)
library(ggplot2)


############################################################
# 1. READ CORTAZAR PTC AND SNV DATA
############################################################

PTC <-
  read.table(
    "/Users/iegab/Downloads/NMD_efficiency_values_PTC.txt",
    header = TRUE
  )


SNV <-
  read.table(
    "/Users/iegab/Downloads/NMD_efficiency_values_SNV.txt",
    header = FALSE
  )


names(SNV) <-
  c(
    "codon",
    "mutation",
    "NMD"
  )


############################################################
# 2. CALCULATE SNV MEDIAN FOR EACH SGE REGION
############################################################

get_region_median <- function(
    data,
    min_codon,
    max_codon
) {
  
  data %>%
    mutate(
      codon = as.numeric(codon),
      NMD = as.numeric(NMD)
    ) %>%
    filter(
      codon >= min_codon,
      codon <= max_codon
    ) %>%
    group_by(
      codon,
      mutation
    ) %>%
    summarise(
      SNVexp_NMD =
        mean(NMD),
      .groups = "drop"
    ) %>%
    summarise(
      SNVmedian_NMD =
        median(SNVexp_NMD)
    ) %>%
    pull(
      SNVmedian_NMD
    )
}


SNV_X1R1 <-
  get_region_median(
    SNV, 1, 27
  )

SNV_X1R2 <-
  get_region_median(
    SNV, 28, 58
  )

SNV_X1R3 <-
  get_region_median(
    SNV, 59, 88
  )

SNV_X1R4 <-
  get_region_median(
    SNV, 89, 118
  )

SNV_X4R1 <-
  get_region_median(
    SNV, 214, 246
  )

SNV_X4R2 <-
  get_region_median(
    SNV, 247, 270
  )

SNV_X7R3 <-
  get_region_median(
    SNV, 434, 460
  )

SNV_X10R1 <-
  get_region_median(
    SNV, 537, 566
  )

SNV_X11R3 <-
  get_region_median(
    SNV, 567, 656
  )


############################################################
# 3. MAP REGIONAL SNV MEDIANS TO PTC CODONS
############################################################

dfx_SNVs_medians <-
  bind_rows(
    
    tibble(
      codon = 1:23,
      SNVmedian_NMD =
        SNV_X1R1
    ),
    
    tibble(
      codon = 24:58,
      SNVmedian_NMD =
        SNV_X1R2
    ),
    
    tibble(
      codon = 59:86,
      SNVmedian_NMD =
        SNV_X1R3
    ),
    
    tibble(
      codon = 87:118,
      SNVmedian_NMD =
        SNV_X1R4
    ),
    
    tibble(
      codon = 214:239,
      SNVmedian_NMD =
        SNV_X4R1
    ),
    
    tibble(
      codon = 240:270,
      SNVmedian_NMD =
        SNV_X4R2
    ),
    
    tibble(
      codon = 434:460,
      SNVmedian_NMD =
        SNV_X7R3
    ),
    
    tibble(
      codon = 537:566,
      SNVmedian_NMD =
        SNV_X10R1
    ),
    
    tibble(
      codon = 567:656,
      SNVmedian_NMD =
        SNV_X11R3
    )
  )


############################################################
# 4. NORMALIZE EXPERIMENTAL PTC NMD
############################################################

PTC_norm <-
  PTC %>%
  select(
    codon,
    mutation,
    mean_val,
    std_val
  ) %>%
  mutate(
    codon =
      as.numeric(codon),
    
    mean_val =
      as.numeric(mean_val),
    
    std_val =
      as.numeric(std_val)
  ) %>%
  
  # Cortazar experimental QC
  filter(
    std_val <= 0.2
  ) %>%
  
  left_join(
    dfx_SNVs_medians,
    by = "codon"
  ) %>%
  
  mutate(
    norm_mean =
      mean_val -
      SNVmedian_NMD
  )


dim(PTC_norm)
# Expected: 763 x 6


summary(
  PTC_norm$norm_mean
)


############################################################
# 5. LOAD 86 SNV-COMPATIBLE STOP-GAIN VARIANTS
############################################################

lmna_snv_final <-
  readRDS(
    "/Users/iegab/Downloads/lmna_snv_final.rds"
  )


############################################################
# 6. ADD NORMALIZED EXPERIMENTAL NMD
############################################################

lmna_validation <-
  lmna_snv_final %>%
  left_join(
    PTC_norm %>%
      select(
        codon,
        stop_codon = mutation,
        SNVmedian_NMD,
        norm_mean
      ),
    by =
      c(
        "codon",
        "stop_codon"
      )
  )


nrow(lmna_validation)
# 86

sum(
  is.na(
    lmna_validation$norm_mean
  )
)

# Expected: 2


############################################################
# 7. IDENTIFY THE TWO QC-FAILED VARIANTS
############################################################

lmna_validation %>%
  filter(
    is.na(norm_mean)
  ) %>%
  select(
    variant_id,
    codon,
    stop_codon,
    experimental_NMD,
    experimental_SD
  )

# Expected:
#
# codon 214 TAG, SD = 0.267
# codon 236 TAG, SD = 0.621


############################################################
# 8. READ TRUNKITTEN PREDICTIONS
############################################################

pred <-
  read.delim(
    "/Users/iegab/Downloads/LMNA_TrunKitten_predictions.tsv"
  )


############################################################
# 9. MERGE PREDICTIONS
############################################################

lmna_validation <-
  lmna_validation %>%
  left_join(
    pred %>%
      select(
        variant_id,
        escape_prob,
        escape_pred_at_youden,
        threshold_used
      ),
    by = "variant_id"
  )


nrow(lmna_validation)
# 86

sum(
  is.na(
    lmna_validation$escape_prob
  )
)

# Expected: 0


############################################################
# 10. FINAL EXPERIMENTAL VALIDATION DATASET
############################################################

lmna_validation_final <-
  lmna_validation %>%
  filter(
    !is.na(norm_mean)
  ) %>%
  mutate(
    
    # Convert escape probability to NMD-oriented probability
    predicted_NMD =
      1 - escape_prob
  )


nrow(lmna_validation_final)
# Expected: 84


############################################################
# 11. PRIMARY SPEARMAN CORRELATION
############################################################

cor.test(
  lmna_validation_final$predicted_NMD,
  lmna_validation_final$norm_mean,
  method = "spearman",
  exact = FALSE
)

# Expected:
#
# rho = 0.5069765
# P = 8.606e-07
# n = 84


############################################################
# 12. MAIN VALIDATION FIGURE
############################################################

ggplot(
  lmna_validation_final,
  aes(
    x = predicted_NMD,
    y = norm_mean,
    fill = factor(mut.exon)
  )
) +
  
  geom_point(
    shape = 21,
    color = "black",
    size = 2.7,
    stroke = 0.5
  ) +
  
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  labs(
    x =
      "TrunKitten-predicted NMD probability",
    
    y =
      "Experimental NMD efficiency",
    
    fill =
      "LMNA exon"
  )


############################################################
# 13. SAVE FINAL AUDIT DATASET
############################################################

write.csv(
  lmna_validation_final,
  "/Users/iegab/Downloads/LMNA_SGE_validation_final.csv",
  row.names = FALSE
)

