library(aenmd)
vcf_file <- "clinvar_20250923.vcf.gz"
vcf <- aenmd:::parse_vcf_VariantAnnotation(vcf_file)
vcf_rng <- vcf$vcf_rng
vcf_rng_fil <- process_variants(vcf_rng)
#- filter out variants with ill-defined alternative allele
ind_out <-  Biostrings::vcountPattern("N", vcf_rng_fil$alt) > 0
vcf_rng_fil <- vcf_rng_fil[!ind_out]
#- back to the original workflow
vcf_rng_ann <- annotate_nmd(vcf_rng_fil, rettype="gr")

saveRDS(vcf_rng_ann, file = "clinvar_aenmd_annotated.rds")

library(S4Vectors)

df <- as.data.frame(vcf_rng_ann, row.names = NULL, optional = TRUE)

# Drop columns that cannot be written to CSV
df <- df[vapply(df, function(x) is.atomic(x) || is.character(x), logical(1))]

