
df$ALLELE.RAT<- -log((df$ALT_COUNT+1)/(df$REF_COUNT+1),2)
df$ALLELE.RAT <- -log((df$altCount+1)/(df$refCount+1),2)

#df$ALLELE.RAT <- (df$REF_COUNT+1)/(df$ALT_COUNT+1)
df$ALLELE.RAT <- (df$REF_COUNT)/(df$ALT_COUNT+df$REF_COUNT)

df$relativePTClocation <- df$coding.pos/df$cds_length

# PTC-bearing exon (rename mut.exon for clarity)
df$PTCBearingExon <- as.integer(df$mut.exon)

# Number of coding exons before and after the PTC-bearing exon
df$AmountExonsBefore <- df$PTCBearingExon - 1
df$AmountExonsAfter  <- df$exon_count - df$PTCBearingExon


out <- mapply(function(pos, exon_idx, cds_exons, exon_len) {

  if (is.na(pos) || is.na(exon_idx) || is.na(cds_exons) || is.na(exon_len)) {
    return(c(NA, NA, NA, NA))
  }

  ends <- as.numeric(strsplit(cds_exons, ",")[[1]])

  if (exon_idx < 1 || exon_idx > length(ends)) {
    return(c(NA, NA, NA, NA))
  }

  exon_end   <- ends[exon_idx]
  exon_start <- if (exon_idx == 1) 1 else ends[exon_idx - 1] + 1

  dist_start_0b <- pos - exon_start
  dist_end_0b   <- exon_end - pos
  dist_start_1b <- dist_start_0b + 1
  rel_exon      <- dist_start_1b / exon_len

  c(dist_start_0b, dist_end_0b, dist_start_1b, rel_exon)

},
pos       = as.numeric(df$coding.pos),
exon_idx  = as.integer(df$mut.exon),
cds_exons = df$cds_exons,
exon_len  = df$length.mutated.exon
)

out <- t(out)
colnames(out) <- c(
  "PTC_dist_exon_start_0b",
  "PTC_dist_exon_end_0b",
  "PTC_dist_exon_start_1b",
  "relPTC_exon"
)

df <- cbind(df, out)

