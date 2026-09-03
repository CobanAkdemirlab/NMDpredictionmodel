# gnomAD_extraction.R
# Title: gnomAD v4.1 extraction
#
# Dataset download/initial parsing follows this script from a collaborator
# repo (confirm which path is current -- two different paths have been used
# in different places):
#   https://github.com/CobanAkdemirlab/NMDescapediseasegene_paper/blob/main/variant%20level_v4/gnomad/gnomAD_downloaddata.R
#
# This script picks up from that repo's output (per-chromosome .rds files)
# and performs: merge -> PTC filtering (via aenmd) -> SNV split -> ANNOVAR.
#
# Requires: config.R, annovar_helpers.R (in this same directory)

source("config.R")
source("annovar_helpers.R")

library(GenomicRanges)
library(aenmd)

# --- 1. Merge per-chromosome gnomAD .rds files (chr1-chrX) ------------------
files <- list.files(cfg$gnomad$raw_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0) {
  stop("No .rds files found in ", cfg$gnomad$raw_dir,
       " -- run gnomAD_downloaddata.R (see header comment) first.")
}

gr_list <- lapply(files, readRDS)
gr_all <- do.call(c, gr_list)
message("Total gnomAD variants merged: ", length(gr_all))

saveRDS(gr_all, out_path(cfg$gnomad, "gnomAD_all_aenmd.rds"))

# --- 2. Identify PTC variants ------------------------------------------------
# FIX: this whole block was missing in the original script -- `df_gnomad_ptc`
# was used below without ever being created. Added here following the same
# process_variants() -> annotate_nmd() pattern used in ClinVar/GREGoR scripts.
# Confirm process_variants() is appropriate to call directly on gr_all here
# (vs. requiring VCF-specific parsing first, as in the other two scripts).
gr_fil <- process_variants(gr_all)
gr_ann <- annotate_nmd(gr_fil, rettype = "gr")

df_gnomad_all <- as.data.frame(gr_ann, row.names = NULL, optional = TRUE)
df_gnomad_ptc <- df_gnomad_all[df_gnomad_all$res_aenmd.is_ptc, ]
message("PTC variants: ", nrow(df_gnomad_ptc))

saveRDS(df_gnomad_ptc, out_path(cfg$gnomad, "gnomAD_isptc.rds"))

# --- 3. Split to SNVs and prepare VCF for ANNOVAR ----------------------------
df_gnomad_snv <- df_gnomad_ptc[df_gnomad_ptc$type == "snv", ]

df_gnomad_snv$key2 <- with(df_gnomad_snv,
  paste0("chr", seqnames, ":", start, "_", ref, ">", alt)
)
df_gnomad_snv$variantID <- with(df_gnomad_snv,
  paste0("chr", seqnames, "_", start, "_", ref, "_", alt)
)

vcf_out_snv <- out_path(cfg$gnomad, "gnomAD_v1_stopgain.vcf")
write_minimal_vcf(df_gnomad_snv, vcf_out_snv)

run_annovar(
  vcf_path    = vcf_out_snv,
  out_prefix  = out_path(cfg$gnomad, paste0("gnomAD_v1_stopgain_gencode_", cfg$gencode_version)),
  annovar_dir = cfg$annovar_dir,
  build       = cfg$genome_build,
  dbtype      = "ensGene"
)

# --- 4. Save canonical-matching table ---------------------------------------
fr.var.can <- data.frame(
  contig     = paste0("chr", df_gnomad_snv$seqnames),
  position   = df_gnomad_snv$start,
  REF_ALLELE = as.character(df_gnomad_snv$ref),
  ALT_ALLELE = as.character(df_gnomad_snv$alt),
  key        = as.character(df_gnomad_snv$key2),
  variantID  = df_gnomad_snv$variantID,
  stringsAsFactors = FALSE
)
save(fr.var.can, file = out_path(cfg$gnomad, "gnomAD_fr.var.can_snv.RData"))

# --- 5. UNRESOLVED: merged_df -----------------------------------------------
# The original script saved an object called `merged_df` at this point:
#   saveRDS(merged_df, file = "merged_chr1to22_snv_pass_df.rds")
#   write.table(merged_df, file = "merged_chr1to22_snv_pass_df.tsv", ...)
# `merged_df` was never defined anywhere in the script you shared, and its
# name/content ("chr1to22", "pass") doesn't obviously match df_gnomad_snv
# (which is chr1-chrX and doesn't reference a PASS filter). This looks like
# it may belong to a different, undocumented step -- e.g. a QC/PASS-filter
# join against something else. Flagging rather than guessing:
#   -> What produces merged_df, and does it feed into the feature matrix,
#      or was it exploratory/unused?

# --- 6. Frameshift variants --------------------------------------------------
# FIX: `df_true` was used in the original script without being defined here.
# Following the pattern from GREGoR/ClinVar, this is the frameshift-type
# subset of the PTC-filtered data. Confirm "frameshift" is the correct
# value in the `type` column for your aenmd version.
df_gnomad_frameshift <- df_gnomad_ptc[df_gnomad_ptc$type == "frameshift", ]

vcf_out_fs <- out_path(cfg$gnomad, "gnomAD_v1_stopgain_frameshift.vcf")
write_minimal_vcf(df_gnomad_frameshift, vcf_out_fs)

run_annovar(
  vcf_path    = vcf_out_fs,
  out_prefix  = out_path(cfg$gnomad, "gnomAD_v1_stopgain_frameshift_gencode_v38"),
  annovar_dir = cfg$annovar_dir,
  build       = cfg$genome_build,
  dbtype      = "ensGene"
)

# Title: Genome Aggregation Database version 4.1 (gnomAD v4.1) extraction
## Dataset extraction followed https://github.com/CobanAkdemirlab/NMDescapediseasegene_paper/blob/main/NMDesc/gnomAD/gnomAD_downloaddata.R
## Next, merge the datasets (chr1-chrX) to select stopgain variants

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


