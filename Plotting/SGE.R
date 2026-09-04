#External validation (SGE)

library(rtracklayer)
library(dplyr)

gtf <- import("gencode.v26.primary_assembly.annotation.gtf")
cds_gtf <- gtf[gtf$type == "CDS"]
cds_df <- as.data.frame(cds_gtf)
cds_df <- cds_df %>%
  mutate(
    CDS_exon_length = end - start + 1
  )
cds_length <- cds_df %>%
  group_by(transcript_id) %>%
  summarise(
    gene_name = dplyr::first(gene_name),
    gene_id = dplyr::first(gene_id),
    strand = dplyr::first(strand),
    cds_length = sum(CDS_exon_length),
    n_coding_exons = dplyr::n(),
    .groups = "drop"
  )
head(cds_length)
cds_length %>%
  filter(gene_name == "LMNA")
