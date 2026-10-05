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
# Paths can either be edited below or supplied on the command line:
#
# Rscript 02_PTBP1_binding_features.R \
#   <GENCODE_GTF> \
#   <OUTPUT_RDS> \
#   <PTBP1_BED1> [<PTBP1_BED2> ...]

args <- commandArgs(
    trailingOnly = TRUE
)

CONFIG <- list(

    gencode_gtf =
        "/path/to/gencode.v26.primary_assembly.annotation.gtf.gz",

    ptbp1_bed_files = c(
        "/path/to/ENCFF907HNN.bed",
        "/path/to/ENCFF130PWU.bed",
        "/path/to/ENCFF100OEX.bed"
    ),

    output_file =
    if (length(args) >= 2) {
        aegs [2]
    } else {
        "/path/to/PTBP1_3UTR_features.rds",
    },
    
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

# Define the proximal 3'UTR as the first 400 nt of the
# spliced 3'UTR in transcript 5' -> 3' orientation.
#
# For transcripts with a 3'UTR shorter than 400 nt, the full
# annotated 3'UTR is retained and is not extended beyond its
# annotated boundary.
#
# For spliced 3'UTRs, the 400 nt are counted cumulatively across
# exonic UTR segments so that intronic sequence is not included.

first_n_of_spliced <- function(
    grl,
    n
) {

    gr <- unlist(
        grl,
        use.names = TRUE
    )


    gr$tx <- names(
        gr
    )


    minus <-
        as.character(
            strand(
                gr
            )
        ) == "-"


    # Order each transcript's UTR segments in transcript
    # 5' -> 3' orientation.
    ord <- order(
        gr$tx,
        ifelse(
            minus,
            -start(gr),
            start(gr)
        )
    )


    gr <- gr[
        ord
    ]


    w <- width(
        gr
    )


    # Number of UTR nucleotides preceding each segment
    # within the same transcript.
    before <- ave(
        w,
        gr$tx,
        FUN = function(x) {
            cumsum(x) - x
        }
    )


    # Number of nucleotides to retain from each segment.
    keep <- pmax(
        0L,
        pmin(
            w,
            n - before
        )
    )


    gr <- gr[
        keep > 0
    ]


    keep <- keep[
        keep > 0
    ]


    # Keep the transcript-oriented start of each UTR segment.
    gr <- resize(
        gr,
        width = keep,
        fix = "start"
    )


    tx <- gr$tx


    mcols(gr) <- NULL
    names(gr) <- NULL


    split(
        gr,
        factor(
            tx,
            levels =
                unique(
                    names(
                        grl
                    )
                )
        )
    )
}


threeutr_proximal <- first_n_of_spliced(
    threeutr_gr,
    CONFIG$proximal_3utr_width
)

# QC:
# each proximal region should contain exactly
# min(400, total annotated 3'UTR length) nucleotides.
stopifnot(
    all(
        sum(
            width(
                threeutr_proximal
            )
        ) ==
        pmin(
            CONFIG$proximal_3utr_width,
            sum(
                width(
                    threeutr_gr
                )
            )
        )
    )
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
    # PTBP1 eCLIP peaks are strand-specific.
    # BED column 6 is therefore required.
    
    if (ncol(bed) < 6) {
        stop(
            "BED file must contain at least three columns: ",
            "including strand in column 6: ",
            bed_file
        )
    }

    peak_strand <- bed[[6]]
    
    peak_strand[
        !peak_strand %in%
            c(
                "+",
                "-"
            )
    ] <- "*"

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
message(
    "PTBP1 peaks by strand: ",
    paste(
        names(
            table(
                strand(
                    ptbp1_gr
                )
            )
        ),
        table(
            strand(
                ptbp1_gr
            )
        ),
        sep = "=",
        collapse = ", "
    )
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

write.table(
    PTBP1_features,
    sub(
        "\\.rds$",
        ".tsv",
        CONFIG$output_file
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
message(
    "PTBP1 3'UTR feature extraction complete."
)

message(
    "Output RDS: ",
    CONFIG$output_file
)


message(
    "Output TSV: ",
    sub(
        "\\.rds$",
        ".tsv",
        CONFIG$output_file
    )
)
