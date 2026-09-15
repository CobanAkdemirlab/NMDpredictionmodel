# ==============================================================================
# PTBP1 3'UTR binding feature
#
# Purpose:
#   Identify transcripts with PTBP1 binding sites overlapping the proximal
#   3'UTR region.
#
# Reference resources:
#   - Genome assembly: GRCh38 / hg38
#   - Transcript annotation: GENCODE v26 primary assembly
#
# External PTBP1 datasets:
#   - ENCFF907HNN
#   - ENCFF130PWU
#   - ENCFF100OEX
#
# Input:
#   - gencode.v26.primary_assembly.annotation.gtf.gz
#   - PTBP1 BED files
#
# Output:
#   - PTBP1_3UTR_features.rds
#
# Output columns:
#   - txnames
#   - threeUTR.PTBP1
#
# This feature table is transcript-level and can subsequently be merged
# with TOPMed, gnomAD, ClinVar, and GREGoR variants by `txnames`.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------------------------

library(GenomicFeatures)
library(GenomicRanges)


# ------------------------------------------------------------------------------
# 2. Configuration
# ------------------------------------------------------------------------------

CONFIG <- list(

    gencode_gtf =
        "/path/to/gencode.v26.primary_assembly.annotation.gtf.gz",

    ptbp1_bed_files = c(
        "/path/to/ENCFF907HNN.bed",
        "/path/to/ENCFF130PWU.bed",
        "/path/to/ENCFF100OEX.bed"
    ),

    output_file =
        "/path/to/PTBP1_3UTR_features.rds",

    # Region of the 3'UTR evaluated in the original analysis
    proximal_3utr_width = 400,

    chromosomes =
        paste0("chr", 1:22)
)


# ------------------------------------------------------------------------------
# 3. Build GENCODE v26 transcript database
# ------------------------------------------------------------------------------

txdb <- makeTxDbFromGFF(
    CONFIG$gencode_gtf
)

txdb <- keepSeqlevels(
    txdb,
    CONFIG$chromosomes,
    pruning.mode = "coarse"
)


# ------------------------------------------------------------------------------
# 4. Retrieve 3'UTR regions
# ------------------------------------------------------------------------------

threeutr_gr <- threeUTRsByTranscript(
    txdb,
    use.names = TRUE
)

message(
    "Transcripts with annotated 3'UTRs: ",
    length(threeutr_gr)
)


# ------------------------------------------------------------------------------
# 5. Define proximal 3'UTR region
# ------------------------------------------------------------------------------

# Reproduce the original analysis by examining the first 400 nt
# of the annotated 3'UTR.
#
# resize(..., fix = "start") respects the strand orientation of the
# GRanges object and extends from the transcript-oriented start.

threeutr_proximal <- resize(
    threeutr_gr,
    width = CONFIG$proximal_3utr_width,
    fix = "start"
)


# ------------------------------------------------------------------------------
# 6. Read PTBP1 binding-site BED files
# ------------------------------------------------------------------------------

read_bed_as_granges <- function(
    bed_file
) {

    bed <- read.table(
        bed_file,
        sep = "\t",
        header = FALSE,
        stringsAsFactors = FALSE
    )

    if (ncol(bed) < 3) {
        stop(
            "BED file must contain at least three columns: ",
            bed_file
        )
    }

    # BED coordinates are 0-based, half-open.
    # GRanges coordinates are 1-based, closed.
    GRanges(
        seqnames = bed[[1]],
        ranges = IRanges(
            start = bed[[2]] + 1,
            end = bed[[3]]
        )
    )
}


ptbp1_list <- lapply(
    CONFIG$ptbp1_bed_files,
    read_bed_as_granges
)

ptbp1_gr <- do.call(
    c,
    ptbp1_list
)

message(
    "PTBP1 binding intervals loaded: ",
    length(ptbp1_gr)
)


# ------------------------------------------------------------------------------
# 7. Harmonize chromosome naming
# ------------------------------------------------------------------------------

common_seqlevels <- intersect(
    seqlevels(threeutr_proximal),
    seqlevels(ptbp1_gr)
)

threeutr_proximal <- keepSeqlevels(
    threeutr_proximal,
    common_seqlevels,
    pruning.mode = "coarse"
)

ptbp1_gr <- keepSeqlevels(
    ptbp1_gr,
    common_seqlevels,
    pruning.mode = "coarse"
)


# ------------------------------------------------------------------------------
# 8. Identify transcripts overlapping PTBP1 binding sites
# ------------------------------------------------------------------------------

ptbp1_hits <- findOverlaps(
    threeutr_proximal,
    ptbp1_gr,
    ignore.strand = TRUE
)

hit_transcripts <- unique(
    names(
        threeutr_proximal
    )[
        queryHits(
            ptbp1_hits
        )
    ]
)

message(
    "Transcripts with proximal 3'UTR PTBP1 overlap: ",
    length(hit_transcripts)
)


# ------------------------------------------------------------------------------
# 9. Construct transcript-level PTBP1 feature
# ------------------------------------------------------------------------------

PTBP1_features <- data.frame(

    txnames =
        names(
            threeutr_gr
        ),

    threeUTR.PTBP1 =
        names(
            threeutr_gr
        ) %in%
        hit_transcripts,

    stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# 10. QC summary
# ------------------------------------------------------------------------------

print(
    table(
        PTBP1_features$threeUTR.PTBP1,
        useNA = "ifany"
    )
)


# ------------------------------------------------------------------------------
# 11. Save output
# ------------------------------------------------------------------------------

saveRDS(
    PTBP1_features,
    CONFIG$output_file
)

message(
    "PTBP1 3'UTR feature extraction complete."
)

message(
    "Output: ",
    CONFIG$output_file
)
