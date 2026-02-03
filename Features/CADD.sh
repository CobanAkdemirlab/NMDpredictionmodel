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

