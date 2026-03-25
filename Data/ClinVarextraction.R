setwd('/home/iegab/TOPMed2026')
library(aenmd)
vcf_file <- "clinvar_20260201.vcf.gz"
vcf <- aenmd:::parse_vcf_VariantAnnotation(vcf_file)
vcf_rng <- vcf$vcf_rng
vcf_rng_fil <- process_variants(vcf_rng)

ind_out <-  Biostrings::vcountPattern("N", vcf_rng_fil$alt) > 0
vcf_rng_fil <- vcf_rng_fil[!ind_out]
vcf_rng_ann <- annotate_nmd(vcf_rng_fil, rettype="gr")

df <- as.data.frame(vcf_rng_ann, row.names = NULL, optional = TRUE)
df_true <- df[df$res_aenmd.is_ptc, ]
df_snv <- df_true[df_true$type == 'snv', ]
save(df_snv, file='/home/iegab/TOPMed2026/clinvar_snv_20260201.RData')
#Create matching columns
df_snv$CHROM <- df_snv$seqnames
df_snv$POS <- df_snv$start
df_snv$key <- with(df_snv,
  paste0("chr", seqnames, ":", start, "_", ref, ">", alt)
)
df_snv$variantID <- with(df_snv,
  paste0("chr", seqnames, "_", start, "_", ref, "_", alt)
)

#prepare the vcf for annovar
vcf <- data.frame(
CHROM = fr.var.can$CHROM,
POS = fr.var.can$start,
ID = ".",
REF = fr.var.can$ref,
ALT = fr.var.can$alt,
QUAL = ".",
FILTER= ".",
INFO = "."
)
#save
write.table(vcf,
file = "/home/iegab/TOPMed2026/clinVar_v1_stopgain.vcf",
sep = "\t", quote = FALSE,row.names = FALSE, col.names = FALSE)
#running annovar
system(paste('perl /home/iegab/annovar/convert2annovar.pl -format vcf4 /home/iegab/TOPMed2026/clinVar_v1_stopgain.vcf > /home/iegab/TOPMed2026/clinVar_v1_stopgain.avinput',sep = ''))
system(paste('perl /home/iegab/annovar/annotate_variation.pl -build hg38 -out /home/iegab/TOPMed2026/clinvar_v1_stopgain_gencode_v38 -dbtype ensGene /home/iegab/TOPMed2026/clinVar_v1_stopgain.avinput /home/iegab/annovar/tempdir', sep = ''))



