# ==============================================================================
# GENCODE v26 sequence feature generation
#
# Purpose:
#   Generate transcript-level sequence-context features used by
#   TrunCat and TrunKitten.
#
# This script generates features from:
#   - coding sequence (CDS)
#   - 5' UTR
#   - 3' UTR
#
# Reference resources:
#   Genome assembly: GRCh38 / hg38
#   Transcript annotation: GENCODE v26 primary assembly
#
# Input:
#   - gencode.v26.primary_assembly.annotation.gtf.gz
#   - BSgenome.Hsapiens.UCSC.hg38
#
# Output:
#   - gencode_v26_sequence_features.rds
#
# Output is keyed by:
#   txnames
#
# This transcript-level feature table is generated once and can be merged
# with annotated TOPMed, gnomAD, ClinVar, and GREGoR variants by txnames.
#
# NOTE:
#   Transcript structural variables required for PTC annotation
#   (cds_length, cds_exons, exon_count, etc.) are generated separately by:
#
#   ../Annotation/02_GENCODEv26_transcript_structure.R
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------------------------

library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(GenomicRanges)


# ------------------------------------------------------------------------------
# 2. Configuration
# ------------------------------------------------------------------------------

CONFIG <- list(

    gencode_gtf =
        "/path/to/gencode.v26.primary_assembly.annotation.gtf.gz",

    output_file =
        "/path/to/gencode_v26_sequence_features.rds",

    # The original analysis used autosomes chr1-chr22.
    chromosomes =
        paste0("chr", 1:22)
)


# ------------------------------------------------------------------------------
# 3. Load GENCODE v26 transcript annotation
# ------------------------------------------------------------------------------

genome <- BSgenome.Hsapiens.UCSC.hg38

txdb <- makeTxDbFromGFF(
    CONFIG$gencode_gtf
)

txdb <- keepSeqlevels(
    txdb,
    CONFIG$chromosomes,
    pruning.mode = "coarse"
)


# ------------------------------------------------------------------------------
# 4. Extract transcript regions
# ------------------------------------------------------------------------------

cds_gr <- cdsBy(
    txdb,
    by = "tx",
    use.names = TRUE
)

threeutr_gr <- threeUTRsByTranscript(
    txdb,
    use.names = TRUE
)

fiveutr_gr <- fiveUTRsByTranscript(
    txdb,
    use.names = TRUE
)

introns_gr <- intronsByTranscript(
    txdb,
    use.names = TRUE
)


# ------------------------------------------------------------------------------
# 5. Extract CDS sequences
# ------------------------------------------------------------------------------

cds_seqs <- extractTranscriptSeqs(
    genome,
    cds_gr
)


# ------------------------------------------------------------------------------
# 6. Extract and collapse UTR sequences by transcript
# ------------------------------------------------------------------------------

collapse_transcript_sequence <- function(
    granges_list,
    genome
) {

    parts <- getSeq(
        genome,
        granges_list
    )

    collapsed <- unlist(
        endoapply(
            parts,
            function(x) {
                do.call(
                    xscat,
                    as.list(x)
                )
            }
        )
    )

    names(collapsed) <-
        names(granges_list)

    collapsed
}


three_utr_seqs <- collapse_transcript_sequence(
    threeutr_gr,
    genome
)

five_utr_seqs <- collapse_transcript_sequence(
    fiveutr_gr,
    genome
)


# ==============================================================================
# Helper functions
# ==============================================================================


# ------------------------------------------------------------------------------
# 7. Sequence length
# ------------------------------------------------------------------------------

get_sequence_length <- function(seq) {

    length(seq)
}


# ------------------------------------------------------------------------------
# 8. Whole-region nucleotide composition
#
# DNA sequence uses T rather than U.
# Therefore:
#   AU content = A + T
#   UC content = T + C
# ------------------------------------------------------------------------------

get_composition <- function(
    seq,
    bases
) {

    freq <- alphabetFrequency(
        seq,
        baseOnly = TRUE,
        as.prob = TRUE
    )

    sum(
        freq[bases]
    )
}


get_AU_content <- function(seq) {

    get_composition(
        seq,
        c("A", "T")
    )
}


get_GC_content <- function(seq) {

    get_composition(
        seq,
        c("G", "C")
    )
}


get_UC_content <- function(seq) {

    get_composition(
        seq,
        c("T", "C")
    )
}


# ------------------------------------------------------------------------------
# 9. First/last sequence-window composition
# ------------------------------------------------------------------------------

get_window_composition <- function(
    seq,
    bases,
    window_size,
    position = c(
        "first",
        "last"
    )
) {

    position <- match.arg(
        position
    )

    if (
        length(seq) <
        window_size
    ) {
        return(
            NA_real_
        )
    }

    if (
        position ==
        "first"
    ) {

        seq_window <- subseq(
            seq,
            start = 1,
            end = window_size
        )

    } else {

        seq_window <- subseq(
            seq,
            start =
                length(seq) -
                window_size +
                1,
            end = length(seq)
        )
    }

    get_composition(
        seq_window,
        bases
    )
}


# ------------------------------------------------------------------------------
# 10. 5'UTR uORF indicator
#
# This reproduces the feature definition used in the original workflow:
# a 5'UTR is considered to contain a candidate uORF when an ATG occurs
# upstream of an in-frame stop codon (TAA, TAG, or TGA).
# ------------------------------------------------------------------------------

has_uORF <- function(seq) {

    seq_string <- as.character(
        seq
    )

    aug_positions <- gregexpr(
        "ATG",
        seq_string
    )[[1]]

    stop_positions <- gregexpr(
        "TAG|TAA|TGA",
        seq_string
    )[[1]]

    aug_positions <-
        aug_positions[
            aug_positions > 0
        ]

    stop_positions <-
        stop_positions[
            stop_positions > 0
        ]

    if (
        length(aug_positions) == 0 ||
        length(stop_positions) == 0
    ) {
        return(
            "not any utr_uORF"
        )
    }

    valid_pair <- any(
        vapply(
            aug_positions,
            function(aug) {

                any(
                    stop_positions > aug &
                    (
                        stop_positions -
                        aug -
                        3
                    ) %% 3 == 0
                )
            },
            logical(1)
        )
    )

    if (valid_pair) {

        "there is a utr_uORF"

    } else {

        "not any utr_uORF"
    }
}


# ------------------------------------------------------------------------------
# 11. Downstream in-frame AUG positions in the CDS
#
# Output is retained under the original feature name:
#   downstream_start
# ------------------------------------------------------------------------------

find_downstream_inframe_aug <- function(seq) {

    seq_string <- as.character(
        seq
    )

    aug_positions <- gregexpr(
        "ATG",
        seq_string
    )[[1]]

    aug_positions <-
        aug_positions[
            aug_positions > 0
        ]

    if (
        length(aug_positions) == 0
    ) {
        return(
            NA_character_
        )
    }

    inframe_aug <- aug_positions[
        aug_positions %% 3 == 1
    ]

    if (
        length(inframe_aug) == 0
    ) {
        return(
            NA_character_
        )
    }

    paste(
        inframe_aug,
        collapse = ","
    )
}


# ------------------------------------------------------------------------------
# 12. Calculate composition features for a DNAStringSet
# ------------------------------------------------------------------------------

calculate_sequence_features <- function(
    seqs,
    prefix
) {

    feature_df <- data.frame(
        txnames = names(seqs),
        stringsAsFactors = FALSE
    )

    # Whole sequence
    feature_df[[paste0(
        prefix,
        "_AU_content"
    )]] <- vapply(
        seqs,
        get_AU_content,
        numeric(1)
    )

    feature_df[[paste0(
        prefix,
        "_GC_content"
    )]] <- vapply(
        seqs,
        get_GC_content,
        numeric(1)
    )

    feature_df[[paste0(
        prefix,
        "_UC_content"
    )]] <- vapply(
        seqs,
        get_UC_content,
        numeric(1)
    )

    # First/last 100 and 200 nucleotides
    for (
        window_size in
        c(100, 200)
    ) {

        for (
            position in
            c("first", "last")
        ) {

            feature_df[[
                paste0(
                    prefix,
                    "_AU_",
                    position,
                    window_size
                )
            ]] <- vapply(
                seqs,
                get_window_composition,
                numeric(1),
                bases =
                    c("A", "T"),
                window_size =
                    window_size,
                position =
                    position
            )

            feature_df[[
                paste0(
                    prefix,
                    "_GC_",
                    position,
                    window_size
                )
            ]] <- vapply(
                seqs,
                get_window_composition,
                numeric(1),
                bases =
                    c("G", "C"),
                window_size =
                    window_size,
                position =
                    position
            )

            feature_df[[
                paste0(
                    prefix,
                    "_UC_",
                    position,
                    window_size
                )
            ]] <- vapply(
                seqs,
                get_window_composition,
                numeric(1),
                bases =
                    c("T", "C"),
                window_size =
                    window_size,
                position =
                    position
            )
        }
    }

    feature_df
}


# ==============================================================================
# CDS features
# ==============================================================================


# ------------------------------------------------------------------------------
# 13. CDS nucleotide composition
# ------------------------------------------------------------------------------

cds_features <- calculate_sequence_features(
    cds_seqs,
    prefix = "CDS"
)


# ------------------------------------------------------------------------------
# 14. Restore original model feature names
#
# The model was trained using these column names, so they are retained
# for compatibility with the original feature matrix.
# ------------------------------------------------------------------------------

names(cds_features)[
    names(cds_features) ==
        "CDS_AU_content"
] <- "cdsseqs_AU_content"

names(cds_features)[
    names(cds_features) ==
        "CDS_GC_content"
] <- "cdsseqs_GC_content"

names(cds_features)[
    names(cds_features) ==
        "CDS_UC_content"
] <- "cdsseqs_UC_content"


# Window names
rename_cds_window <- c(

    CDS_AU_last100 =
        "cdsseq_AUcontentlast100",

    CDS_GC_last100 =
        "cdsseq_GCcontentlast100",

    CDS_UC_last100 =
        "cdsseq_UCcontentlast100",

    CDS_AU_first100 =
        "cdsseq_AUcontentfirst100",

    CDS_GC_first100 =
        "cdsseq_GCcontentfirst100",

    CDS_UC_first100 =
        "cdsseq_UCcontentfirst100",

    CDS_AU_last200 =
        "cdsseq_AUcontentlast200",

    CDS_GC_last200 =
        "cdsseq_GCcontentlast200",

    CDS_UC_last200 =
        "cdsseq_UCcontentlast200",

    CDS_AU_first200 =
        "cdsseq_AUcontentfirst200",

    CDS_GC_first200 =
        "cdsseq_GCcontentfirst200",

    CDS_UC_first200 =
        "cdsseq_UCcontentfirst200"
)

for (
    old_name in
    names(rename_cds_window)
) {

    names(cds_features)[
        names(cds_features) ==
            old_name
    ] <-
        rename_cds_window[
            old_name
        ]
}


# ------------------------------------------------------------------------------
# 15. Downstream in-frame AUG positions
# ------------------------------------------------------------------------------

cds_features$downstream_start <- vapply(
    cds_seqs[
        cds_features$txnames
    ],
    find_downstream_inframe_aug,
    character(1)
)


# ==============================================================================
# 3'UTR features
# ==============================================================================


# ------------------------------------------------------------------------------
# 16. 3'UTR length and sequence composition
# ------------------------------------------------------------------------------

threeutr_features <-
    calculate_sequence_features(
        three_utr_seqs,
        prefix = "ThreeUTR"
    )

threeutr_features$threeUTR_length <-
    vapply(
        three_utr_seqs[
            threeutr_features$txnames
        ],
        get_sequence_length,
        integer(1)
    )


# ------------------------------------------------------------------------------
# 17. Restore original 3'UTR feature names
# ------------------------------------------------------------------------------

names(threeutr_features)[
    names(threeutr_features) ==
        "ThreeUTR_AU_content"
] <- "threeUTR_AU_content"

names(threeutr_features)[
    names(threeutr_features) ==
        "ThreeUTR_GC_content"
] <- "threeUTR_GC_content"

names(threeutr_features)[
    names(threeutr_features) ==
        "ThreeUTR_UC_content"
] <- "threeUTR_UC_content"


rename_threeutr_window <- c(

    ThreeUTR_AU_last100 =
        "ThreeUTR_AUcontentlast100",

    ThreeUTR_GC_last100 =
        "ThreeUTR_GCcontentlast100",

    ThreeUTR_UC_last100 =
        "ThreeUTR_UCcontentlast100",

    ThreeUTR_AU_first100 =
        "ThreeUTR_AUcontentfirst100",

    ThreeUTR_GC_first100 =
        "ThreeUTR_GCcontentfirst100",

    ThreeUTR_UC_first100 =
        "ThreeUTR_UCcontentfirst100",

    ThreeUTR_AU_last200 =
        "ThreeUTR_AUcontentlast200",

    ThreeUTR_GC_last200 =
        "ThreeUTR_GCcontentlast200",

    ThreeUTR_UC_last200 =
        "ThreeUTR_UCcontentlast200",

    ThreeUTR_AU_first200 =
        "ThreeUTR_AUcontentfirst200",

    ThreeUTR_GC_first200 =
        "ThreeUTR_GCcontentfirst200",

    ThreeUTR_UC_first200 =
        "ThreeUTR_UCcontentfirst200"
)

for (
    old_name in
    names(rename_threeutr_window)
) {

    names(threeutr_features)[
        names(threeutr_features) ==
            old_name
    ] <-
        rename_threeutr_window[
            old_name
        ]
}


# ==============================================================================
# 5'UTR features
# ==============================================================================


# ------------------------------------------------------------------------------
# 18. 5'UTR length and sequence composition
# ------------------------------------------------------------------------------

fiveutr_features <-
    calculate_sequence_features(
        five_utr_seqs,
        prefix = "FiveUTR"
    )

fiveutr_features$fiveutr_length <-
    vapply(
        five_utr_seqs[
            fiveutr_features$txnames
        ],
        get_sequence_length,
        integer(1)
    )


# ------------------------------------------------------------------------------
# 19. 5'UTR uORF feature
# ------------------------------------------------------------------------------

fiveutr_features$fiveutrseqs.uORF <-
    vapply(
        five_utr_seqs[
            fiveutr_features$txnames
        ],
        has_uORF,
        character(1)
    )


# ------------------------------------------------------------------------------
# 20. Restore original 5'UTR feature names
# ------------------------------------------------------------------------------

names(fiveutr_features)[
    names(fiveutr_features) ==
        "FiveUTR_AU_content"
] <- "fiveUTR_AU_content"

names(fiveutr_features)[
    names(fiveutr_features) ==
        "FiveUTR_GC_content"
] <- "fiveUTR_GC_content"

names(fiveutr_features)[
    names(fiveutr_features) ==
        "FiveUTR_UC_content"
] <- "fiveUTR_UC_content"


rename_fiveutr_window <- c(

    FiveUTR_AU_last100 =
        "fiveUTR_AUcontentlast100",

    FiveUTR_GC_last100 =
        "fiveUTR_GCcontentlast100",

    FiveUTR_UC_last100 =
        "fiveUTR_UCcontentlast100",

    FiveUTR_AU_first100 =
        "fiveUTR_AUcontentfirst100",

    FiveUTR_GC_first100 =
        "fiveUTR_GCcontentfirst100",

    FiveUTR_UC_first100 =
        "fiveUTR_UCcontentfirst100",

    FiveUTR_AU_last200 =
        "FiveUTR_AUcontentlast200",

    FiveUTR_GC_last200 =
        "FiveUTR_GCcontentlast200",

    FiveUTR_UC_last200 =
        "FiveUTR_UCcontentlast200",

    FiveUTR_AU_first200 =
        "FiveUTR_AUcontentfirst200",

    FiveUTR_GC_first200 =
        "FiveUTR_GCcontentfirst200",

    FiveUTR_UC_first200 =
        "FiveUTR_UCcontentfirst200"
)

for (
    old_name in
    names(rename_fiveutr_window)
) {

    names(fiveutr_features)[
        names(fiveutr_features) ==
            old_name
    ] <-
        rename_fiveutr_window[
            old_name
        ]
}


# ==============================================================================
# UTR intron indicators
# ==============================================================================


# ------------------------------------------------------------------------------
# 21. Transcript-level UTR intron status
#
# The feature indicates whether the transcript contains more than one
# segment of the corresponding UTR, consistent with an intron interrupting
# that UTR.
# ------------------------------------------------------------------------------

has_multiple_regions <- function(
    gr_list
) {

    vapply(
        gr_list,
        length,
        integer(1)
    ) > 1
}


threeutr_intron_table <- data.frame(

    txnames =
        names(
            threeutr_gr
        ),

    threeUTR.introns =
        ifelse(
            has_multiple_regions(
                threeutr_gr
            ),
            "There is a 3UTR intron",
            "NA"
        ),

    stringsAsFactors = FALSE
)


fiveutr_intron_table <- data.frame(

    txnames =
        names(
            fiveutr_gr
        ),

    fiveUTR.introns =
        ifelse(
            has_multiple_regions(
                fiveutr_gr
            ),
            "There is a 5UTR intron",
            "NA"
        ),

    stringsAsFactors = FALSE
)


# ==============================================================================
# Merge transcript-level sequence features
# ==============================================================================


# ------------------------------------------------------------------------------
# 22. Merge CDS, 3'UTR, and 5'UTR features
# ------------------------------------------------------------------------------

gencode_v26_sequence_features <-
    merge(
        cds_features,
        threeutr_features,
        by = "txnames",
        all.x = TRUE
    )

gencode_v26_sequence_features <-
    merge(
        gencode_v26_sequence_features,
        fiveutr_features,
        by = "txnames",
        all.x = TRUE
    )

gencode_v26_sequence_features <-
    merge(
        gencode_v26_sequence_features,
        threeutr_intron_table,
        by = "txnames",
        all.x = TRUE
    )

gencode_v26_sequence_features <-
    merge(
        gencode_v26_sequence_features,
        fiveutr_intron_table,
        by = "txnames",
        all.x = TRUE
    )


# ------------------------------------------------------------------------------
# 23. Save transcript-level feature table
# ------------------------------------------------------------------------------

saveRDS(
    gencode_v26_sequence_features,
    CONFIG$output_file
)

message(
    "GENCODE v26 sequence feature extraction complete."
)

message(
    "Transcripts: ",
    nrow(
        gencode_v26_sequence_features
    )
)

message(
    "Feature columns: ",
    ncol(
        gencode_v26_sequence_features
    ) - 1
)

  
  
