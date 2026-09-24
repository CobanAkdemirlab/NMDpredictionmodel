############################################################
# 02_LMNA_add_conservation_and_predict.R
#
# Purpose:
#   Merge conservation features with the prepared LMNA
#   variants and construct the exact 10-feature input
#   required by TrunKitten.
#
# Inputs:
#   1. lmna_snv_before_conservation.rds
#   2. variants_with_conservation_medians.csv
#
# Output:
#   LMNA_TrunKitten_10features_with_id.tsv
############################################################
suppressPackageStartupMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste0(
      "Usage:\n",
      "Rscript 02_LMNA_add_conservation_and_predict.R ",
      "<LMNA_RDS> <CONSERVATION_CSV> <OUTPUT_TSV>"
    )
  )
}

LMNA_FILE <- args[1]
CONSERVATION_FILE <- args[2]
OUTPUT_FILE <- args[3]

lmna_snv <- readRDS(LMNA_FILE)
conservation <- read.csv(CONSERVATION_FILE)

############################################################
# ADD CONSERVATION FEATURES
############################################################

library(dplyr)

# Reload intermediate object if needed
lmna_snv <-
  readRDS(
    "/path/lmna_snv_before_conservation.rds"
  )


conservation <-
  read.csv(
    "/path/variants_with_conservation_medians.csv"
  )


lmna_snv_final <-
  lmna_snv %>%
  left_join(
    conservation %>%
      select(
        variant_id,
        phastcons_new3utr_first200_median,
        phylop_ptc_to_ejc_median
      ),
    by = "variant_id"
  )


# QC
nrow(lmna_snv_final)

sum(
  is.na(
    lmna_snv_final$
      phastcons_new3utr_first200_median
  )
)

sum(
  is.na(
    lmna_snv_final$
      phylop_ptc_to_ejc_median
  )
)

# Expected:
#
# 86
# 0
# 0


############################################################
# CREATE EXACT 10-FEATURE TRUNKITTEN INPUT
############################################################

trunkitten_input_with_id <-
  lmna_snv_final %>%
  select(
    variant_id,
    last.EJC,
    relativePTClocation,
    half_life_PC1,
    cdsseqs_AU_content,
    mut.exon,
    phastcons_new3utr_first200_median,
    phylop_ptc_to_ejc_median,
    AmountExonsAfter,
    cdsseq_AUcontentlast200,
    cdsseqs_UC_content
  )


dim(trunkitten_input_with_id)

# Expected:
# 86 rows x 11 columns
# 1 ID + 10 features


write.table(
  trunkitten_input_with_id,
  "LMNA_TrunKitten_10features_with_id.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message(
  "TrunKitten input written: ",
  OUTPUT_FILE
)

#then on ssh
cd /home/iegab/TrunKitten_LMNA

source trunkitten_env/bin/activate
python NMDpredictionmodel/Model/TrunKitten/pipeline/predict.py \
--annotated /home/iegab/TrunKitten_LMNA/LMNA_TrunKitten_10features_with_id.tsv \
--model /home/iegab/TrunKitten_LMNA/NMDpredictionmodel/Model/TrunKitten/model/trunkitten.pkl \
--metadata /home/iegab/TrunKitten_LMNA/NMDpredictionmodel/Model/TrunKitten/model/trunkitten_features.json \
--out /home/iegab/TrunKitten_LMNA/LMNA_TrunKitten_predictions.tsv

#scored 86 variants
#predicted escape: 45 (52.3%)
#threshold: 0.4218

#proceed to 03_LMNA_TrunKitten_validation.R
