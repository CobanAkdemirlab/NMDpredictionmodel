library(GenomicRanges)
gr_all <- do.call(c, rds.merged)
gr_all
length(gr_all)
names(gr_all) <- make.unique(names(gr_all))
gregor <- as.data.frame(gr_all)
#6127
#prepare for ANNOVAR
vcf_gregor <- data.frame(
    CHROM  = paste0("chr", gregor$seqnames),
    POS    = gregor$start,
    ID     = ".",
    REF    = as.character(gregor$ref),
    ALT    = as.character(gregor$alt),
    QUAL   = ".",
    FILTER = "PASS",
    INFO   = ".",
    stringsAsFactors = FALSE
)
