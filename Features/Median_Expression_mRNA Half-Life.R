#ENSG 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
library(biomaRt)
mart <- useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
BM.info <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","hgnc_symbol","transcript_is_canonical"),mart=mart)
BM.info.can <- BM.info[which(BM.info$transcript_is_canonical=='1'),]
#df$TxName <- sub("\\..*$","", df$txnames)

df <- merge(df, BM.info.can, by.x='txnames', by.y='ensembl_transcript_id', all.x=TRUE)


# Median Expression
exp <- read.table('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.txt',header=TRUE,sep='\t')
exp$gene_id <- sub("\\..*", "", exp$Name)

# identify tissue columns (everything except Name, Description, gene_id)
tissue_cols <- setdiff(colnames(exp), c("Name", "Description", "gene_id"))

# compute median expression per gene
exp$MedianExpression <- apply(
  exp[, tissue_cols],
  1,
  median,
  na.rm = TRUE
)

df.merged <- merge(
  df,
  exp[, c("gene_id", "MedianExpression", "Whole.Blood")],
  by.x = "ensembl_gene_id",
  by.y = "gene_id",
  all.x = TRUE
)

df$MedianExpression_log2 <- log2(df$MedianExpression + 1)

# mRNAHalfLifeMin (PC)
library(readxl)


file <- "/location/13059_2022_2811_MOESM3_ESM.xlsx"

# list sheets first
excel_sheets(file)
# skip first 2 rows so row 3 becomes the header
hl3 <- read_excel(file, sheet = "human", skip = 1)


#merging
df_merged <- merge(
    df,
    hl3[, c("Ensembl Gene Id", "half-life (PC1)")],
    by.x = "ensembl_gene_id",
    by.y = "Ensembl Gene Id",
    all.x = TRUE
)
#fix the `half-life (PC1)` if necessary
df$half_life_PC1 <- df$`half-life (PC1)`
#scaling half-life PC1 based on TOPMed
mean_ref <- mean(df$half_life_PC1, na.rm = TRUE) #TOPMed (used as a reference)
sd_ref   <- sd(df$half_life_PC1, na.rm = TRUE) #TOPMed (used as a reference)
#then apply to datasets
df$half_life_pc1_z <- (df$half_life_PC1 - mean_ref) / sd_ref #change this based in the propor df TOPMed, gnomAD, clinVar and GREGoR

breaks <- quantile(df$half_life_pc1_z,
                   probs = c(0, 1/3, 2/3, 1),
                   na.rm = TRUE)
df_merged <- df_merged %>%
    mutate(
        HL_bin = cut(
            half_life_pc1_z,
            breaks = breaks,
            labels = c("short", "medium", "long"),
            include.lowest = TRUE
        )
    )
