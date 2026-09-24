############################################################
# 01_LMNA_prepare_TrunKitten_features.R
#
# Purpose:
#   Prepare SNV-compatible stop-gain variants from the
#   Cortazar et al. LMNA saturation genome editing dataset
#   and reconstruct the non-conservation TrunKitten features.
#
# Inputs:
#   1. Cortazar et al. LMNA PTC measurements
#   2. GENCODE v26 GTF
#   3. GRCh38 reference FASTA
#   4. mRNA half-life PC1 table
#
# Outputs:
#   LMNA_Cortazar_TrunKitten_variants.csv
#   lmna_snv_before_conservation.rds
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(rtracklayer)
  library(GenomicRanges)
  library(Rsamtools)
  library(Biostrings)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop(
    paste0(
      "Usage:\n",
      "Rscript 01_LMNA_prepare_TrunKitten_features.R ",
      "<PTC_FILE> <GENCODE_V26_GTF> <GRCh38_FASTA> ",
      "<HALF_LIFE_PC1_XLSX> <OUTPUT_DIR>"
    )
  )
}

PTC_FILE      <- args[1]
GTF_FILE      <- args[2]
GENOME_FILE   <- args[3]
HALF_LIFE_FILE <- args[4]
OUTPUT_DIR    <- args[5]

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

LMNA_TX <- "ENST00000368300.8"

############################################################
# 1. INPUT FILES
############################################################

PTC_FILE <- "/path/NMD_efficiency_values_PTC.txt"

GTF_FILE <- "/path/gencode.v26.primary_assembly.annotation.gtf"

GENOME_FILE <- "/path/GRCh38.primary_assembly.genome.fa"

HALF_LIFE_FILE <- "/path/half_life_pc1.xlsx"

# Transcript used for LMNA
LMNA_TX <- "ENST00000368300.8"


############################################################
# 2. READ CORTAZAR PTC DATA
############################################################

PTC <- read.table(
  PTC_FILE,
  header = TRUE
)

dim(PTC)
head(PTC)


############################################################
# 3. CREATE INITIAL LMNA DATASET
############################################################

lmna <- PTC %>%
  transmute(
    codon,
    stop_codon = mutation,
    experimental_NMD = mean_val,
    experimental_SD = std_val,
    gene = "LMNA"
  )

nrow(lmna)
# Expected: 770


############################################################
# 4. READ GENCODE v26 AND SELECT LMNA TRANSCRIPT CDS
############################################################

gtf <- import(GTF_FILE)

cds_gtf <- gtf[gtf$type == "CDS"]

cds_df <- as.data.frame(cds_gtf)

lmna_tx <- cds_df %>%
  filter(
    gene_name == "LMNA",
    transcript_id == LMNA_TX
  ) %>%
  arrange(start) %>%
  mutate(
    CDS_exon_length = end - start + 1,
    CDS_cum_end = cumsum(CDS_exon_length),
    CDS_cum_start = CDS_cum_end - CDS_exon_length + 1
  )

lmna_tx %>%
  select(
    seqnames,
    start,
    end,
    strand,
    transcript_id,
    exon_number,
    CDS_exon_length,
    CDS_cum_start,
    CDS_cum_end
  )

LMNA_CDS_LENGTH <- sum(lmna_tx$CDS_exon_length)

LMNA_CDS_LENGTH
# Expected: 1992


############################################################
# 5. CONVERT CODON NUMBER -> CDS CODING POSITION
############################################################

# Position of first nucleotide of each codon:
#
# codon 1 -> 1
# codon 2 -> 4
# codon 3 -> 7
# etc.

lmna <- lmna %>%
  mutate(
    coding.pos = 3 * (codon - 1) + 1
  )


############################################################
# 6. CALCULATE relativePTClocation
############################################################

lmna <- lmna %>%
  mutate(
    relativePTClocation =
      coding.pos / LMNA_CDS_LENGTH
  )


############################################################
# 7. ASSIGN mut.exon
############################################################

get_coding_exon <- function(coding_pos, tx) {
  
  hit <- tx %>%
    filter(
      coding_pos >= CDS_cum_start,
      coding_pos <= CDS_cum_end
    )
  
  if (nrow(hit) == 0) {
    return(NA_integer_)
  }
  
  as.integer(hit$exon_number[1])
}


lmna <- lmna %>%
  rowwise() %>%
  mutate(
    mut.exon = get_coding_exon(
      coding.pos,
      lmna_tx
    )
  ) %>%
  ungroup()


table(lmna$mut.exon)

# Expected experimental exons:
# 1, 4, 7, 10, 11


############################################################
# 8. AmountExonsAfter
############################################################

# LMNA transcript contains 12 transcript exons.

TOTAL_LMNA_EXONS <- 12

lmna <- lmna %>%
  mutate(
    AmountExonsAfter =
      TOTAL_LMNA_EXONS - mut.exon
  )


lmna %>%
  distinct(
    mut.exon,
    AmountExonsAfter
  ) %>%
  arrange(mut.exon)

# Expected:
#
# exon 1  -> 11
# exon 4  -> 8
# exon 7  -> 5
# exon 10 -> 2
# exon 11 -> 1


############################################################
# 9. LAST EXON RULE
############################################################

lmna <- lmna %>%
  mutate(
    last.exon = if_else(
      mut.exon == TOTAL_LMNA_EXONS,
      "lastexon",
      "notlastexon"
    )
  )


############################################################
# 10. PENULTIMATE LAST-50-BP RULE
############################################################

# Reproduce the original TrunCat definition.
#
# We use cumulative CDS exon ends.

cds_exon_ends <- lmna_tx$CDS_cum_end


lmna$penultimate.last50bp <- sapply(
  1:nrow(lmna),
  function(x) {
    
    coding_pos <- lmna$coding.pos[x]
    mut_exon <- lmna$mut.exon[x]
    
    exon_count <- TOTAL_LMNA_EXONS
    
    # Check whether variant is in penultimate exon
    if (exon_count == mut_exon + 1) {
      
      penultimate_end <-
        cds_exon_ends[length(cds_exon_ends) - 1]
      
      penultimate_start <-
        cds_exon_ends[length(cds_exon_ends) - 2]
      
      penultimate_length <-
        penultimate_end - penultimate_start
      
      fifty_pos <-
        penultimate_end - 50
      
      if (
        coding_pos >= fifty_pos &&
        coding_pos <= penultimate_end &&
        penultimate_length >= 50
      ) {
        
        return("penultimate.last50bp")
      }
    }
    
    return("not penultimate.last50bp")
  }
)


############################################################
# 11. CREATE last.EJC
############################################################

lmna <- lmna %>%
  mutate(
    last.EJC = case_when(
      last.exon == "lastexon" ~
        "last.exon",
      
      penultimate.last50bp ==
        "penultimate.last50bp" ~
        "penultimate.last50bp",
      
      TRUE ~
        "upstream"
    )
  )


table(lmna$last.EJC)

# For Cortazar LMNA:
#
# expected only:
# upstream
# penultimate.last50bp
#
# No last exon was experimentally assayed.


############################################################
# 12. READ LMNA CDS SEQUENCE FROM GRCh38
############################################################

genome <- FaFile(GENOME_FILE)

open(genome)

lmna_cds_gr <- GRanges(
  seqnames = lmna_tx$seqnames,
  ranges = IRanges(
    start = lmna_tx$start,
    end = lmna_tx$end
  ),
  strand = lmna_tx$strand
)

cds_parts <- getSeq(
  genome,
  lmna_cds_gr
)

LMNA_CDS_SEQ <-
  paste0(
    as.character(cds_parts),
    collapse = ""
  )

LMNA_CDS_SEQ <- toupper(LMNA_CDS_SEQ)

nchar(LMNA_CDS_SEQ)
# Expected: 1992


############################################################
# 13. SEQUENCE COMPOSITION FEATURES
############################################################

seq <- LMNA_CDS_SEQ

L <- nchar(seq)

A_count <- str_count(seq, "A")
T_count <- str_count(seq, "T")
C_count <- str_count(seq, "C")
G_count <- str_count(seq, "G")


LMNA_AU_content <-
  (A_count + T_count) / L

LMNA_UC_content <-
  (T_count + C_count) / L


LMNA_CDS_last200 <-
  substr(
    seq,
    L - 199,
    L
  )

LMNA_AU_last200 <-
  (
    str_count(LMNA_CDS_last200, "A") +
      str_count(LMNA_CDS_last200, "T")
  ) / 200


LMNA_AU_content
# ~0.3664659

LMNA_UC_content
# ~0.4442771

LMNA_AU_last200
# ~0.34


lmna <- lmna %>%
  mutate(
    cdsseqs_AU_content =
      LMNA_AU_content,
    
    cdsseqs_UC_content =
      LMNA_UC_content,
    
    cdsseq_AUcontentlast200 =
      LMNA_AU_last200
  )


############################################################
# 14. ADD LMNA HALF-LIFE PC1
############################################################

half_life <-
  read_excel(HALF_LIFE_FILE)


LMNA_half_life_PC1 <-
  half_life %>%
  filter(`Gene name` == "LMNA") %>%
  pull(half_life_PC1)


LMNA_half_life_PC1
# Expected ~4.190347


lmna <- lmna %>%
  mutate(
    half_life_PC1 =
      LMNA_half_life_PC1
  )


############################################################
# 15. RECONSTRUCT NATIVE LMNA CODON
############################################################

lmna <- lmna %>%
  mutate(
    ref_codon = str_sub(
      LMNA_CDS_SEQ,
      coding.pos,
      coding.pos + 2
    )
  )


lmna %>%
  distinct(
    codon,
    coding.pos,
    ref_codon
  ) %>%
  arrange(codon) %>%
  head(15)


############################################################
# 16. COUNT NUMBER OF NUCLEOTIDE CHANGES
############################################################

count_diffs <- function(ref, alt) {
  
  ref_chars <-
    strsplit(ref, "")[[1]]
  
  alt_chars <-
    strsplit(alt, "")[[1]]
  
  sum(
    ref_chars != alt_chars
  )
}


lmna <- lmna %>%
  rowwise() %>%
  mutate(
    n_nt_changes =
      count_diffs(
        ref_codon,
        stop_codon
      )
  ) %>%
  ungroup()


table(lmna$n_nt_changes)

# Expected:
#
# 1 change = 86
# 2 changes = 288
# 3 changes = 396
#
# Total = 770


############################################################
# 17. KEEP TRUE SINGLE-NUCLEOTIDE STOP-GAIN VARIANTS
############################################################

lmna_snv <- lmna %>%
  filter(
    n_nt_changes == 1
  )


nrow(lmna_snv)
# Expected: 86


############################################################
# 18. IDENTIFY EXACT ALTERED BASE WITHIN CODON
############################################################

find_changed_base <- function(ref, alt) {
  
  which(
    strsplit(ref, "")[[1]] !=
      strsplit(alt, "")[[1]]
  )[1]
}


get_base <- function(codon, pos) {
  
  substr(
    codon,
    pos,
    pos
  )
}


lmna_snv <- lmna_snv %>%
  rowwise() %>%
  mutate(
    changed_base_in_codon =
      find_changed_base(
        ref_codon,
        stop_codon
      ),
    
    variant_coding_pos =
      coding.pos +
      changed_base_in_codon -
      1,
    
    refAllele =
      get_base(
        ref_codon,
        changed_base_in_codon
      ),
    
    altAllele =
      get_base(
        stop_codon,
        changed_base_in_codon
      )
  ) %>%
  ungroup()


############################################################
# 19. CONVERT CDS POSITION -> GENOMIC POSITION
############################################################

coding_to_genomic <- function(cdna_pos) {
  
  hit <- lmna_tx %>%
    filter(
      cdna_pos >= CDS_cum_start,
      cdna_pos <= CDS_cum_end
    )
  
  if (nrow(hit) == 0) {
    return(NA_real_)
  }
  
  hit$start[1] +
    (
      cdna_pos -
        hit$CDS_cum_start[1]
    )
}


lmna_snv <- lmna_snv %>%
  rowwise() %>%
  mutate(
    position =
      coding_to_genomic(
        variant_coding_pos
      )
  ) %>%
  ungroup()


############################################################
# 20. CREATE STANDARD VARIANT IDENTIFIERS
############################################################

lmna_snv <- lmna_snv %>%
  mutate(
    contig = "chr1",
    gene = "LMNA",
    txnames = LMNA_TX,
    strand = "+",
    
    variant_id =
      paste0(
        contig, ":",
        position, ":",
        refAllele, ">",
        altAllele
      )
  )


############################################################
# 21. QC: CONFIRM REF ALLELES AGAINST GRCh38
############################################################

test_gr <- GRanges(
  seqnames = lmna_snv$contig,
  ranges = IRanges(
    start = lmna_snv$position,
    end = lmna_snv$position
  )
)

genome_ref <-
  as.character(
    getSeq(
      genome,
      test_gr
    )
  )


table(
  expected = lmna_snv$refAllele,
  genome = genome_ref
)

# Expected: perfect diagonal match
# for all 86 variants.


############################################################
# 22. EXPORT VARIANTS FOR CONSERVATION ANNOTATION
############################################################

LMNA_Cortazar_TrunKitten_variants <-
  lmna_snv %>%
  select(
    variant_id,
    contig,
    position,
    refAllele,
    altAllele,
    gene,
    txnames,
    strand
  )


write.csv(
  LMNA_Cortazar_TrunKitten_variants,
  file.path(
    OUTPUT_DIR,
    "LMNA_Cortazar_TrunKitten_variants.csv"
  ),
  row.names = FALSE
)

saveRDS(
  lmna_snv,
  file.path(
    OUTPUT_DIR,
    "lmna_snv_before_conservation.rds"
  )
)

close(genome)

message("LMNA TrunKitten feature preparation complete.")
