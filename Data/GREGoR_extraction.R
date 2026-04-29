#aenmd
library(aenmd)
library(GenomicRanges)
library(GenomeInfoDb)


vcf_file <- "filtered_oc_base_matches.unique.vcf"
vcf <- aenmd:::parse_vcf_VariantAnnotation(vcf_file)
vcf_rng <- vcf$vcf_rng
seqlevels(vcf_rng) <- sub("^chr", "", seqlevels(vcf_rng))
seqnames(vcf_rng) <- sub("^chr", "", as.character(seqnames(vcf_rng)))
vcf_rng_fil <- process_variants(vcf_rng)
length(vcf_rng_fil)

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
write.table(
    vcf_gregor,
    file = "/home/iegab/TOPMed2026/GREGoR_v1_stopgain.vcf",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
)
system(paste(
    'perl ~/annovar/convert2annovar.pl -format vcf4',
    '~/GREGoR_v1_stopgain_frameshift.vcf',
    '>',
    '~/GREGoR_v1_stopgain_frameshift.avinput' ,sep = ''
))

system(paste(
    'perl /home/iegab/annovar/annotate_variation.pl',
    '-build hg38',
    '-out /home/iegab/TOPMed2026/GREGoR_v1_stopgain_hg38',
    '-dbtype ensGene',
    '/home/iegab/TOPMed2026/GREGoR_v1_stopgain.avinput',
    '/home/iegab/annovar/tempdir', sep = ''
))


