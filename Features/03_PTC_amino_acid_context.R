# ==============================================================================
# PTC amino-acid context features
#
# Purpose:
#   Generate amino-acid sequence context features surrounding each
#   premature termination codon (PTC).
#
# Features:
#   - PTC.minus2.aa : amino acid two residues upstream of the PTC
#   - PTC.minus1.aa : amino acid immediately upstream of the PTC
#   - PTC.plus1.aa  : amino acid immediately downstream of the PTC
#   - PTC.plus2.aa  : amino acid two residues downstream of the PTC
#   - PTC.id        : resulting stop-codon identity for stop-gain SNVs
#
# Reference resources:
#   - Genome assembly: GRCh38 / hg38
#   - Transcript annotation: GENCODE v26 primary assembly
#
# Input:
#   - annotated PTC SNV dataset
#   - gencode.v26.primary_assembly.annotation.gtf.gz
#
# Output:
#   - PTC_amino_acid_features.rds
#
# This workflow is shared across TOPMed, gnomAD, ClinVar, and GREGoR.
#
# NOTE:
#   The final study workflow is restricted to SNVs. Frameshift-specific
#   processing from the exploratory analysis is therefore not included here.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------------------------

library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(dplyr)
library(stringr)


# ------------------------------------------------------------------------------
# 2. Configuration
# ------------------------------------------------------------------------------

CONFIG <- list(

    gencode_gtf =
        "/path/to/gencode.v26.primary_assembly.annotation.gtf.gz",

    input_file =
        "/path/to/PTC_annotated_variants.rds",

    output_file =
        "/path/to/PTC_amino_acid_features.rds",

    chromosomes =
        paste0("chr", 1:22)
)


# ------------------------------------------------------------------------------
# 3. Load annotated variants
# ------------------------------------------------------------------------------

variants <- readRDS(
    CONFIG$input_file
)


# ------------------------------------------------------------------------------
# 4. Build GENCODE v26 transcript database
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
# 5. Extract and translate CDS sequences
# ------------------------------------------------------------------------------

genome <- BSgenome.Hsapiens.UCSC.hg38

cds_by_tx <- cdsBy(
    txdb,
    by = "tx",
    use.names = TRUE
)

cds_seqs <- extractTranscriptSeqs(
    genome,
    cds_by_tx
)

protein_seqs <- translate(
    cds_seqs
)


# ------------------------------------------------------------------------------
# 6. Extract protein position from ANNOVAR annotation
# ------------------------------------------------------------------------------

extract_protein_position <- function(
    annotation,
    transcript_id
) {

    fields <- strsplit(
        annotation,
        ":|,"
    )[[1]]

    tx_index <- grep(
        transcript_id,
        fields,
        fixed = TRUE
    )

    if (length(tx_index) == 0) {
        return(NA_real_)
    }

    protein_field <- fields[
        tx_index[1] + 3
    ]

    if (
        length(protein_field) == 0 ||
        !grepl("p\\.", protein_field)
    ) {
        return(NA_real_)
    }

    # For complex annotation strings, retain the first coordinate.
    protein_field <- strsplit(
        protein_field,
        "_"
    )[[1]][1]

    protein_position <- suppressWarnings(
        as.numeric(
            gsub(
                "\\D",
                "",
                protein_field
            )
        )
    )

    protein_position
}


variants$protein_position <- mapply(
    extract_protein_position,
    variants$V3,
    variants$txnames
)


# ------------------------------------------------------------------------------
# 7. Retrieve amino acid relative to PTC
# ------------------------------------------------------------------------------

get_relative_amino_acid <- function(
    transcript_id,
    protein_position,
    offset,
    protein_sequences
) {

    if (
        is.na(transcript_id) ||
        is.na(protein_position)
    ) {
        return(NA_character_)
    }

    tx_index <- match(
        transcript_id,
        names(protein_sequences)
    )

    if (is.na(tx_index)) {
        return(NA_character_)
    }

    sequence <- as.character(
        protein_sequences[
            tx_index
        ]
    )

    target_position <-
        protein_position +
        offset

    if (
        target_position < 1 ||
        target_position >
        nchar(sequence)
    ) {
        return(NA_character_)
    }

    substr(
        sequence,
        target_position,
        target_position
    )
}


# ------------------------------------------------------------------------------
# 8. Generate amino-acid context features
# ------------------------------------------------------------------------------

variants$PTC.minus2.aa <- mapply(
    get_relative_amino_acid,
    variants$txnames,
    variants$protein_position,
    MoreArgs = list(
        offset = -2,
        protein_sequences = protein_seqs
    )
)

variants$PTC.minus1.aa <- mapply(
    get_relative_amino_acid,
    variants$txnames,
    variants$protein_position,
    MoreArgs = list(
        offset = -1,
        protein_sequences = protein_seqs
    )
)

variants$PTC.plus1.aa <- mapply(
    get_relative_amino_acid,
    variants$txnames,
    variants$protein_position,
    MoreArgs = list(
        offset = 1,
        protein_sequences = protein_seqs
    )
)

variants$PTC.plus2.aa <- mapply(
    get_relative_amino_acid,
    variants$txnames,
    variants$protein_position,
    MoreArgs = list(
        offset = 2,
        protein_sequences = protein_seqs
    )
)


# ------------------------------------------------------------------------------
# 9. Determine resulting stop-codon identity
# ------------------------------------------------------------------------------

get_ptc_identity <- function(
    transcript_id,
    annotation,
    cds_sequence
) {

    fields <- strsplit(
        annotation,
        ":|,"
    )[[1]]

    tx_index <- grep(
        transcript_id,
        fields,
        fixed = TRUE
    )

    if (length(tx_index) == 0) {
        return(NA_character_)
    }

    coding_field <- fields[
        tx_index[1] + 2
    ]

    if (!grepl("c\\.", coding_field)) {
        return(NA_character_)
    }

    # Capture the changed nucleotide from the coding annotation.
    alt_base <- substr(
        coding_field,
        nchar(coding_field),
        nchar(coding_field)
    )

    coding_position <- suppressWarnings(
        as.numeric(
            gsub(
                "\\D",
                "",
                coding_field
            )
        )
    )

    if (is.na(coding_position)) {
        return(NA_character_)
    }

    if (coding_position %% 3 == 1) {

        codon_start <- coding_position
        changed_position <- 1

    } else if (
        coding_position %% 3 == 2
    ) {

        codon_start <- coding_position - 1
        changed_position <- 2

    } else {

        codon_start <- coding_position - 2
        changed_position <- 3
    }

    codon <- substr(
        cds_sequence,
        codon_start,
        codon_start + 2
    )

    if (nchar(codon) != 3) {
        return(NA_character_)
    }

    substr(
        codon,
        changed_position,
        changed_position
    ) <- alt_base

    codon
}


variants$PTC.id <- mapply(
    function(
        transcript_id,
        annotation
    ) {

        tx_index <- match(
            transcript_id,
            names(cds_seqs)
        )

        if (is.na(tx_index)) {
            return(NA_character_)
        }

        get_ptc_identity(
            transcript_id =
                transcript_id,

            annotation =
                annotation,

            cds_sequence =
                as.character(
                    cds_seqs[
                        tx_index
                    ]
                )
        )
    },

    variants$txnames,
    variants$V3
)


# ------------------------------------------------------------------------------
# 10. QC
# ------------------------------------------------------------------------------

message(
    "Variants processed: ",
    nrow(variants)
)

message(
    "Variants with PTC -1 amino acid: ",
    sum(
        !is.na(
            variants$PTC.minus1.aa
        )
    )
)

message(
    "Stop codon distribution:"
)

print(
    table(
        variants$PTC.id,
        useNA = "ifany"
    )
)


# ------------------------------------------------------------------------------
# 11. Save feature table
# ------------------------------------------------------------------------------

PTC_amino_acid_features <- variants %>%
    select(
        key,
        variantID,
        txnames,
        protein_position,
        PTC.minus2.aa,
        PTC.minus1.aa,
        PTC.plus1.aa,
        PTC.plus2.aa,
        PTC.id
    )

saveRDS(
    PTC_amino_acid_features,
    CONFIG$output_file
)

message(
    "PTC amino-acid feature extraction complete."
)


#categorization
#PTC (amino acid)
df$minus1_is_G <- rep('Not Glycine',nrow(df))
ind.1 <- which(df$PTC.minus1.aa=='G')
df$minus1_is_G[ind.1] <- rep('Glycine',length(ind.1))


df$minus2_is_G <- rep('Not Glycine',nrow(df))
ind.1 <- which(df$PTC.minus2.aa=='G')
df$minus2_is_G[ind.1] <- rep('Glycine',length(ind.1))

df$plus1_is_G <- rep('Not Glycine',nrow(df))
ind.1 <- which(df$PTC.plus1.aa=='G')
df$plus1_is_G[ind.1] <- rep('Glycine',length(ind.1))

df$plus2_is_G <- rep('Not Glycine',nrow(df))
ind.1 <- which(df$PTC.plus2.aa=='G')
df$plus2_is_G[ind.1] <- rep('Glycine',length(ind.1))


