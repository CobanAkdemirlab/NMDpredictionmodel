# ==============================================================================
# TOPMed ASE simulation for recurrent PTC variants
#
# Purpose:
# Construct one representative allele-specific expression (ASE) ratio
# per annotated TOPMed PTC variant for model development.
#
#   Variants observed in one individual retain their observed ASE value.
#   Variants observed in multiple individuals are repeatedly sampled across
#   carriers and summarized using the median simulated ASE value.
#
# Input:
#   1. Annotated TOPMed variant table (`df`)
#   2. Per-individual ASE/genotype table (`fr.var.can`)
#
# Output:
#   `df.sim` - one record per variant with representative ASE and NMD class
#
# Filtering:
#   totalCount >= 8
#
# Simulation:
#   100 iterations for recurrent variants
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------------------------

library(dplyr)


# ------------------------------------------------------------------------------
# 2. Configuration
# ------------------------------------------------------------------------------

CONFIG <- list(
    annotation_file = "/path/to/TOPMed_stopgain_fulldf_feb1826.RData",
    min_total_count = 8,
    n_simulations = 100,
    random_seed = 1234
)


# ------------------------------------------------------------------------------
# 3. Load annotated TOPMed dataset
# ------------------------------------------------------------------------------

load(CONFIG$annotation_file)

# IMPORTANT:
# The loaded file must create the annotated variant object `df`.
# If possible, replace load() with readRDS() in the final workflow so that
# the object name is explicit.


# ------------------------------------------------------------------------------
# 4. Filter ASE observations by sequencing depth
# ------------------------------------------------------------------------------

fr.var.can.filtered <- fr.var.can %>%
    filter(
        totalCount >= CONFIG$min_total_count
    )

message(
    "ASE observations after depth filtering: ",
    nrow(fr.var.can.filtered)
)


# ------------------------------------------------------------------------------
# 5. Restrict annotation and ASE tables to shared variants
# ------------------------------------------------------------------------------

common_variants <- intersect(
    df$key,
    fr.var.can.filtered$key
)

fr.var.can.filtered <- fr.var.can.filtered %>%
    filter(
        key %in% common_variants
    )

df.filtered <- df %>%
    filter(
        key %in% common_variants
    )

message(
    "Annotated variants with matching ASE data: ",
    length(common_variants)
)


# ------------------------------------------------------------------------------
# 6. Identify unique and recurrent variants
# ------------------------------------------------------------------------------

variant_frequency <- fr.var.can.filtered %>%
    count(
        variantID,
        name = "n_carriers"
    )

unique_variants <- variant_frequency %>%
    filter(
        n_carriers == 1
    ) %>%
    pull(
        variantID
    )

recurrent_variants <- variant_frequency %>%
    filter(
        n_carriers > 1
    ) %>%
    pull(
        variantID
    )

message(
    "Variants observed in one individual: ",
    length(unique_variants)
)

message(
    "Variants observed in multiple individuals: ",
    length(recurrent_variants)
)


# ------------------------------------------------------------------------------
# 7. Unique variants: retain observed ASE
# ------------------------------------------------------------------------------

df.unique <- df.filtered %>%
    filter(
        variantID %in% unique_variants
    ) %>%
    left_join(
        fr.var.can.filtered %>%
            filter(
                variantID %in% unique_variants
            ) %>%
            select(
                variantID,
                refCount,
                altCount,
                totalCount
            ),
        by = "variantID"
    ) %>%
    mutate(
        ALLELE.RAT = refCount / totalCount
    )


# ------------------------------------------------------------------------------
# 8. Recurrent variants: Monte Carlo carrier sampling
# ------------------------------------------------------------------------------

set.seed(CONFIG$random_seed)

simulated_iterations <- vector(
    "list",
    CONFIG$n_simulations
)

for (t in seq_len(CONFIG$n_simulations)) {

    message(
        "Simulation ",
        t,
        " of ",
        CONFIG$n_simulations
    )

    simulated_iteration <- df.filtered %>%
        filter(
            variantID %in% recurrent_variants
        )

    simulated_iteration$ALLELE.RAT <- NA_real_

    for (i in seq_len(nrow(simulated_iteration))) {

        carrier_rows <- which(
            fr.var.can.filtered$variantID ==
                simulated_iteration$variantID[i]
        )

        sampled_row <- sample(
            carrier_rows,
            size = 1
        )

        sampled_ref <- fr.var.can.filtered$refCount[
            sampled_row
        ]

        sampled_total <- fr.var.can.filtered$totalCount[
            sampled_row
        ]

        simulated_iteration$ALLELE.RAT[i] <-
            sampled_ref / sampled_total
    }

    simulated_iterations[[t]] <- simulated_iteration
}


# ------------------------------------------------------------------------------
# 9. Summarize recurrent-variant simulations
# ------------------------------------------------------------------------------

df.sim.common <- bind_rows(
    simulated_iterations
) %>%
    group_by(
        key
    ) %>%
    mutate(
        ALLELE.RAT = median(
            ALLELE.RAT,
            na.rm = TRUE
        )
    ) %>%
    distinct(
        key,
        .keep_all = TRUE
    ) %>%
    ungroup()


# ------------------------------------------------------------------------------
# 10. Combine unique and recurrent variants
# ------------------------------------------------------------------------------

df.sim <- bind_rows(
    df.unique,
    df.sim.common
)

# ------------------------------------------------------------------------------
# 11. Define NMD outcome from representative allele ratio
# ------------------------------------------------------------------------------

# Allele ratio definition:
#
#   ALLELE.RAT = refCount / (refCount + altCount)
#
# Since totalCount = refCount + altCount:
#
#   ALLELE.RAT = refCount / totalCount
#
# For variants observed in one individual, ALLELE.RAT is calculated
# directly from the observed read counts.
#
# For recurrent variants, ALLELE.RAT is calculated within each carrier
# sampling iteration and the median across 100 iterations is retained as
# the representative allele ratio.

df.sim <- df.sim %>%
    mutate(
        NMD.ESCAPEE = case_when(

            ALLELE.RAT >= 0.35 &
                ALLELE.RAT <= 0.65 ~ "TRUE",

            ALLELE.RAT > 0.65 ~ "FALSE",

            TRUE ~ NA_character_
        )
    )
# Allele-specific expression outcome:
#   ALLELE.RAT = refCount / (refCount + altCount)
#
# NMD classification:
#   0.35 <= ALLELE.RAT <= 0.65  -> NMD escape
#   ALLELE.RAT > 0.65           -> NMD-sensitive / stronger NMD
#   ALLELE.RAT < 0.35           -> excluded from binary outcome
# ------------------------------------------------------------------------------
# 12. Recalculate NMD class after ASE aggregation
# ------------------------------------------------------------------------------

df.sim <- df.sim %>%
    mutate(
        NMD.ESCAPEE = case_when(
            ALLELE.RAT >= 0.35 &
                ALLELE.RAT <= 0.65 ~ "TRUE",

            ALLELE.RAT > 0.65 ~ "FALSE",

            TRUE ~ NA_character_
        )
    )


# ------------------------------------------------------------------------------
# 13. QC summary
# ------------------------------------------------------------------------------

message(
    "Final simulated dataset variants: ",
    nrow(df.sim)
)

print(
    table(
        df.sim$NMD.ESCAPEE,
        useNA = "ifany"
    )
)
