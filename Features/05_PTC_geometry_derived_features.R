# ==============================================================================
# Derived PTC geometry features
#
# Purpose:
#   Generate transcript-relative and exon-relative PTC position features.
#
# Shared across:
#   TOPMed, gnomAD, ClinVar, GREGoR
# ==============================================================================

library(dplyr)


# ------------------------------------------------------------------------------
# 1. Load annotated variants
# ------------------------------------------------------------------------------

variants <- readRDS(
    "/path/to/PTC_annotated_variants.rds"
)


# ------------------------------------------------------------------------------
# 2. Relative PTC position within CDS
# ------------------------------------------------------------------------------

variants <- variants %>%
    mutate(
        relativePTClocation =
            coding.pos / cds_length,

        PTCBearingExon =
            as.integer(mut.exon),

        AmountExonsBefore =
            PTCBearingExon - 1,

        AmountExonsAfter =
            exon_count -
            PTCBearingExon
    )


# ------------------------------------------------------------------------------
# 3. Exon-relative PTC geometry
# ------------------------------------------------------------------------------

get_exon_relative_geometry <- function(
    pos,
    exon_idx,
    cds_exons,
    exon_len
) {

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

    ends <- as.numeric(
        strsplit(
            cds_exons,
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

    exon_end <- ends[
        exon_idx
    ]

    exon_start <-
        if (
            exon_idx == 1
        ) {
            1
        } else {
            ends[
                exon_idx - 1
            ] + 1
        }

    dist_start_0b <-
        pos -
        exon_start

    dist_end_0b <-
        exon_end -
        pos

    dist_start_1b <-
        dist_start_0b + 1

    rel_exon <-
        dist_start_1b /
        exon_len

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
}


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
        variants$length.mutated.exon
)

geometry <- t(
    geometry
)

variants <- cbind(
    variants,
    geometry
)
