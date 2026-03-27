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
#after 

