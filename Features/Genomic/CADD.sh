#prepare the dataset for annovar
annovar_input <- data.frame(
  chr   = variants_min$chr,
  start = variants_min$pos,
  end   = variants_min$pos,
  ref   = variants_min$ref,
  alt   = variants_min$alt,
  stringsAsFactors = FALSE
)

#save
write.table(
  annovar_input,
  file = "topmed_for_cadd.avinput",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


#on terminal
perl table_annovar.pl \
  topmed_for_cadd.avinput \
  humandb/ \
  -buildver hg38 \
  -out topmed_cadd \
  -remove \
  -protocol dbnsfp42a,gnomad_exome \
  -operation f,f \
  -nastring . \
  -otherinfo
#merge by key (Chr, Start, End, Ref, 	Alt, REVEL_score, CADD_raw, CADD_phred, gnomAD_exome_ALL== chr7:127588544_A>T)
cadd$key <- with(cadd,
  paste0("chr", Chr, ":", Start, "_", Ref, ">", Alt)
)
cadd_to_merge <- cadd[, c(
  "key",
  "REVEL_score",
  "CADD_raw",
  "CADD_phred",
  "gnomAD_exome_ALL"
)]
df.merged <- merge(
  df,
  cadd_to_merge,
  by = "key",
  all.x = TRUE
)

