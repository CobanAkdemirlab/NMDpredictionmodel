# ==============================================================================
# Derived PTC geometry features
#
# Purpose:
#   Generate transcript-relative and exon-relative positional features for
#   premature termination codon (PTC) variants.
#
# Shared across:
#   - TOPMed
#   - gnomAD
#   - ClinVar
#   - GREGoR
#
# Required input columns:
#   - key
#   - coding.pos
#   - cds_length
#   - mut.exon
#   - exon_count
#   - cds_exons
#   - length.mutated.exon
#
# Features generated:
#   - relativePTClocation
#   - PTCBearingExon
#   - AmountExonsBefore
#   - AmountExonsAfter
#   - PTC_dist_exon_start_0b
#   - PTC_dist_exon_end_0b
#   - PTC_dist_exon_start_1b
#   - relPTC_exon
#
# Output:
#   - PTC_geometry_features.rds
#
# The resulting feature table is keyed by `key` and is merged into the final
# model feature matrix in 99_build_feature_matrix.R.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------------------------

library(dplyr)


# ------------------------------------------------------------------------------
# 2. Configuration
# ------------------------------------------------------------------------------

CONFIG <- list(

    input_file =
        "/path/to/PTC_annotated_variants.rds",

    output_file =
        "/path/to/PTC_geometry_features.rds"
)


# ------------------------------------------------------------------------------
# 3. Load annotated variants
# ------------------------------------------------------------------------------

variants <- readRDS(
    CONFIG$input_file
)


# ------------------------------------------------------------------------------
# 4. Check required columns
# ------------------------------------------------------------------------------

required_columns <- c(
    "key",
    "coding.pos",
    "cds_length",
    "mut.exon",
    "exon_count",
    "cds_exons",
    "length.mutated.exon"
)

missing_columns <- setdiff(
    required_columns,
    colnames(variants)
)

if (length(missing_columns) > 0) {

    stop(
        "Missing required columns: ",
        paste(
            missing_columns,
            collapse = ", "
        )
    )
}


# ------------------------------------------------------------------------------
# 5. Transcript-relative PTC geometry
# ------------------------------------------------------------------------------

variants <- variants %>%
    mutate(

        # Relative position of the PTC within the coding sequence.
        #
        # 0 -> near CDS start
        # 1 -> near CDS end
        relativePTClocation =
            as.numeric(coding.pos) /
            as.numeric(cds_length),

        # Coding exon containing the PTC
        PTCBearingExon =
            as.integer(mut.exon),

        # Number of coding exons before the PTC-bearing exon
        AmountExonsBefore =
            PTCBearingExon - 1L,

        # Number of coding exons downstream of the PTC-bearing exon
        AmountExonsAfter =
            as.integer(exon_count) -
            PTCBearingExon
    )


# ------------------------------------------------------------------------------
# 6. Exon-relative PTC geometry
# ------------------------------------------------------------------------------

get_exon_relative_geometry <- function(
    pos,
    exon_idx,
    cds_exons,
    exon_len
) {

    # --------------------------------------------------------------------------
    # Missing input
    # --------------------------------------------------------------------------

    if (
        is.na(pos) ||
        is.na(exon_idx) ||
        is.na(cds_exons) ||
        is.na(exon_len)
    ) {

        return(
            c(
                PTC_dist_exon_start_0b = NA_real_,
                PTC_dist_exon_end_0b   = NA_real_,
                PTC_dist_exon_start_1b = NA_real_,
                relPTC_exon            = NA_real_
            )
        )
    }


    # --------------------------------------------------------------------------
    # Parse cumulative CDS exon-end coordinates
    #
    # Example:
    #
    #   cds_exons = "120,245,390"
    #
    # means:
    #   exon 1: CDS positions   1-120
    #   exon 2: CDS positions 121-245
    #   exon 3: CDS positions 246-390
    # --------------------------------------------------------------------------

    ends <- as.numeric(
        strsplit(
            as.character(cds_exons),
            ","
        )[[1]]
    )


    if (
        exon_idx < 1 ||
        exon_idx > length(ends)
    ) {

        return(
            c(
                PTC_dist_exon_start_0b = NA_real_,
                PTC_dist_exon_end_0b   = NA_real_,
                PTC_dist_exon_start_1b = NA_real_,
                relPTC_exon            = NA_real_
            )
        )
    }


    # --------------------------------------------------------------------------
    # Determine CDS coordinates of the PTC-bearing exon
    # --------------------------------------------------------------------------

    exon_end <- ends[
        exon_idx
    ]

    exon_start <-
        if (exon_idx == 1) {

            1

        } else {

            ends[
                exon_idx - 1
            ] + 1
        }


    # --------------------------------------------------------------------------
    # Distance from beginning of PTC-bearing exon
    #
    # 0-based:
    #   first nucleotide in exon -> 0
    #
    # 1-based:
    #   first nucleotide in exon -> 1
    # --------------------------------------------------------------------------

    dist_start_0b <-
        pos -
        exon_start

    dist_start_1b <-
        dist_start_0b + 1


    # --------------------------------------------------------------------------
    # Distance from PTC to end of PTC-bearing exon
    #
    # A value of 0 means that the PTC falls at the final CDS nucleotide of
    # the exon.
    # --------------------------------------------------------------------------

    dist_end_0b <-
        exon_end -
        pos


    # --------------------------------------------------------------------------
    # Relative PTC position within PTC-bearing exon
    #
    # Values near:
    #   0 -> exon start
    #   1 -> exon end
    # --------------------------------------------------------------------------

    rel_exon <-
        dist_start_1b /
        exon_len


    return(
        c(
            PTC_dist_exon_start_0b =
                dist_start_0b,

            PTC_dist_exon_end_0b =
                dist_end_0b,

            PTC_dist_exon_start_1b =
                dist_start_1b,

            relPTC_exon =
                rel_exon
        )
    )
}


# ------------------------------------------------------------------------------
# 7. Calculate exon-relative geometry for all variants
# ------------------------------------------------------------------------------

geometry <- mapply(

    get_exon_relative_geometry,

    pos =
        as.numeric(
            variants$coding.pos
        ),

    exon_idx =
        as.integer(
            variants$mut.exon
        ),

    cds_exons =
        variants$cds_exons,

    exon_len =
        as.numeric(
            variants$length.mutated.exon
        )
)


geometry <- t(
    geometry
)


# Ensure numeric output
geometry <- as.data.frame(
    geometry,
    stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# 8. Add geometry features to variant table
# ------------------------------------------------------------------------------

variants <- bind_cols(
    variants,
    geometry
)


# ------------------------------------------------------------------------------
# 9. Quality-control checks
# ------------------------------------------------------------------------------

# Relative position within CDS should generally be between 0 and 1.
invalid_relative_cds <- which(
    !is.na(variants$relativePTClocation) &
    (
        variants$relativePTClocation < 0 |
        variants$relativePTClocation > 1
    )
)

if (length(invalid_relative_cds) > 0) {

    warning(
        length(invalid_relative_cds),
        " variants have relativePTClocation outside [0,1]."
    )
}


# Relative position within the PTC-bearing exon should also generally
# be between 0 and 1.
invalid_relative_exon <- which(
    !is.na(variants$relPTC_exon) &
    (
        variants$relPTC_exon < 0 |
        variants$relPTC_exon > 1
    )
)

if (length(invalid_relative_exon) > 0) {

    warning(
        length(invalid_relative_exon),
        " variants have relPTC_exon outside [0,1]."
    )
}


# Downstream exon count should not be negative.
invalid_exon_count <- which(
    !is.na(variants$AmountExonsAfter) &
    variants$AmountExonsAfter < 0
)

if (length(invalid_exon_count) > 0) {

    warning(
        length(invalid_exon_count),
        " variants have negative AmountExonsAfter."
    )
}


# ------------------------------------------------------------------------------
# 10. QC summary
# ------------------------------------------------------------------------------

message(
    "PTC geometry feature generation complete."
)

message(
    "Variants processed: ",
    nrow(variants)
)

message(
    "Variants with relative CDS position: ",
    sum(
        !is.na(
            variants$relativePTClocation
        )
    )
)

message(
    "Variants with exon-relative position: ",
    sum(
        !is.na(
            variants$relPTC_exon
        )
    )
)


# ------------------------------------------------------------------------------
# 11. Create feature-only output table
# ------------------------------------------------------------------------------

PTC_geometry_features <- variants %>%
    select(
        key,

        relativePTClocation,

        PTCBearingExon,

        AmountExonsBefore,

        AmountExonsAfter,

        PTC_dist_exon_start_0b,

        PTC_dist_exon_end_0b,

        PTC_dist_exon_start_1b,

        relPTC_exon
    )


# ------------------------------------------------------------------------------
# 12. Save
# ------------------------------------------------------------------------------

saveRDS(
    PTC_geometry_features,
    CONFIG$output_file
)

message(
    "Saved: ",
    CONFIG$output_file
)
