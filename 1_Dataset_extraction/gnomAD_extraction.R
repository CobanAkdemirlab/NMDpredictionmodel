# Title: Genome Aggregation Database version 4.1 (gnomAD v4.1) extraction
# Description: 
#gnomAD v4.1
library(aenmd)
library(GenomicRanges)
library(GenomeInfoDb)
Run_Annotate_and_Save <- function(vcf_rng_proc, outputname) {
	  vcf_rng_an  <- aenmd::annotate_nmd(vcf_rng_proc,rettype="gr")
  saveRDS(vcf_rng_an, file= paste0(outputname, '.rds'))
    #return(vcf_rng_an)
}
chrs <- c(2:22)
for (chr in chrs) {
	print(chr)
	chr.no <- paste('chr',chr,sep='')
	link=paste0('https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.',chr.no,'.vcf.bgz')
	system(paste0('wget ',link,sep=' '))

	### get the header of a vcf file

	infile <- paste0("gnomad.exomes.v4.1.sites.",chr.no,'.vcf.bgz')
	truncfile <- paste0("gnomad.exomes.v4.1.sites.",chr.no,'.cut.trunc.vcf')
	trunc1file <- paste0("gnomad.exomes.v4.1.sites.",chr.no,'.cut.trunc.1.vcf')

	system(paste0("gunzip -c ", infile, " | grep '#' > header.vcf"))

	system(paste0(
		        "gunzip -c ", infile,
			  " | grep -v '^#'",
			  " | awk -F'\t' '$7==\"PASS\"'",
			    " | grep -E 'stop_gained|frameshift' > ",
			    truncfile
			    ))

	system(paste0("cat header.vcf ", truncfile, " > " , trunc1file))


	vcf_file <- trunc1file
	vcf      <- aenmd:::parse_vcf_VariantAnnotation(vcf_file,param = VariantAnnotation::ScanVcfParam(geno=NA))
	outputname <- chr.no

	vcf_rng <- vcf$vcf_rng
	vcf_rng$ref     <- vcf_rng$ref |> Biostrings::DNAStringSet()
	vcf_rng$alt     <- vcf_rng$alt |> Biostrings::DNAStringSet()

	seqlevels(vcf_rng) <- gsub('chr','',seqlevels(vcf_rng))
	seqnames(vcf_rng) <- gsub('chr','',seqnames(vcf_rng))

	seqlevelsStyle(vcf_rng) <- 'NCBI' #gives you 1.....X
	#seqlevelsStyle(vcf_rng) <- 'UCSC' #gives you chr1. ....X

	#Funciton to make unique IDs
	make_keys <- function(vcf_rng){
		  starts <- GenomicRanges::start(vcf_rng) |> stringr::str_pad(9L, pad="0")
	  keys <- paste0(GenomicRanges::seqnames(vcf_rng), ":", starts,"|" ,vcf_rng$ref, "|", vcf_rng$alt)
	    return(keys)
	}
	#Make keys, which are unique IDs
	vcf_rng$key <- aenmd:::make_keys(vcf_rng)
	#vcf_ifo <- VariantAnnotation::info(vcf)

	genome(vcf_rng) <- gsub('gnomAD_','',genome(vcf_rng))

	vcf_rng_proc <- aenmd::process_variants(vcf_rng)


	Run_Annotate_and_Save(vcf_rng_proc,outputname)
	rm(vcffile)
	file.remove(infile)
	file.remove(truncfile)
	#- split ranges and info
	#vcf_rng <- SummarizedExperiment::rowRanges(vcf)
	#colnames( vcf_rng |> S4Vectors::mcols() ) <- vcf_rng |> S4Vectors::mcols() |> 
	  #colnames() |> janitor::make_clean_names()

	#Some variants have the alt of "<DEL>" for some weird reason, just make those alt alleles empty:
	#vcf_rng$alt <- ifelse(vcf_rng$alt == "<DEL>", "", vcf_rng$alt)
	#make sure that ref and alt are class biostring:




}

setwd("~/Datasets/gnomAD")
files <- list.files(pattern = "\\.rds$", full.names = TRUE)
gr_list <- lapply(files, readRDS)
library(GenomicRanges)

gr_all <- do.call(c, gr_list)
gr_all
length(gr_all)
head(mcols(gr_all))
saveRDS(gr_all, "~/Datasets/gnomAD/gnomAD_all_aenmd.rds")
  
df_gnomad_snv <- df_gnomad_ptc[which(df_gnomad_ptc$type=='snv'),]
vcf_gnomad <- data.frame(
    CHROM  = paste0("chr", df_gnomad_snv$seqnames),
    POS    = df_gnomad_snv$start,
    ID     = ".",
    REF    = as.character(df_gnomad_snv$ref),
    ALT    = as.character(df_gnomad_snv$alt),
    QUAL   = ".",
    FILTER = "PASS",
    INFO   = ".",
    stringsAsFactors = FALSE
)
writeLines(
    c("##fileformat=VCFv4.2",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"),
    "~/gnomAD/gnomAD_v1_stopgain.vcf"
)
write.table(
    vcf_gnomad,
    file = "~/gnomAD/gnomAD_v1_stopgain.vcf",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
)
system(paste(
    'perl ~/annovar/convert2annovar.pl -format vcf4',
    '~/gnomAD/gnomAD_v1_stopgain.vcf',
    '>',
    '~/gnomAD/gnomAD_v1_stopgain.avinput'
))
system(paste(
    'perl ~/annovar/annotate_variation.pl',
    '-build hg38',
    '-out ~/gnomAD/gnomAD_v1_stopgain_gencode_v38',
    '-dbtype ensGene',
    '~/gnomAD/gnomAD_v1_stopgain.avinput',
    '~/annovar/tempdir'
))
fr.var.can <- data.frame(
    contig      = paste0("chr", df_gnomad_snv$seqnames),
    position    = df_gnomad_snv$start,
    REF_ALLELE  = as.character(df_gnomad_snv$ref),
    ALT_ALLELE  = as.character(df_gnomad_snv$alt),
    key = as.character(df_gnomad_snv$key2),
    variantID = df_gnomad_snv$variantID,
    stringsAsFactors = FALSE
)
save(fr.var.can, file = "/home/iegab/TOPMed2026/gnomAD/gnomAD_fr.var.can_snv.RData")
# Optionally save 

saveRDS(merged_df, file = "merged_chr1to22_snv_pass_df.rds") 

write.table(merged_df, file = "merged_chr1to22_snv_pass_df.tsv", sep = "\t", row.names = FALSE, quote = FALSE) 
vcf <- data.frame(
    CHROM = paste0("chr", df_true$seqnames),
    POS   = df_true$start,
    ID    = ".",
    REF   = df_true$ref,
    ALT   = df_true$alt,
    QUAL  = ".",
    FILTER= ".",
    INFO  = "."
)

write.table(vcf,
          file = "~/gnomAD/gnomAD_v1_stopgain_frameshift.vcf",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
system(paste('perl ~/annovar/convert2annovar.pl -format vcf4 ~/gnomAD/gnomAD_v1_stopgain_frameshift.vcf > ~/gnomAD/gnomAD_v1_stopgain_frameshift.avinput',sep = ''))

system(paste('perl ~/annovar/annotate_variation.pl -build hg38 -out ~/gnomAD/gnomAD_v1_stopgain_frameshift_gencode_v38 -dbtype ensGene ~/gnomAD/gnomAD_v1_stopgain_frameshift.avinput ~/annovar/tempdir', sep = ''))


