# ==============================================================================
# Prepare variants for external score annotation
#
# Purpose:
#   Convert PTC SNVs into ANNOVAR input format for annotation with
#   dbNSFP/CADD and gnomAD population allele frequencies.
#
# Shared across:
#   TOPMed, gnomAD, ClinVar, GREGoR
#
# Output:
#   variants_for_external_annotation.avinput
# ==============================================================================

CONFIG <- list(
    input_file =
        "/path/to/annotated_variants.rds",

    output_file =
        "/path/to/variants_for_external_annotation.avinput"
)

variants <- readRDS(
    CONFIG$input_file
)

# ------------------------------------------------------------------------------
# Prepare ANNOVAR input
# ------------------------------------------------------------------------------

annovar_input <- data.frame(
    chr   = variants$chr,
    start = variants$pos,
    end   = variants$pos,
    ref   = variants$ref,
    alt   = variants$alt,
    stringsAsFactors = FALSE
)

write.table(
    annovar_input,
    file = CONFIG$output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

message(
    "ANNOVAR input written: ",
    CONFIG$output_file
)

message(
    "Variants: ",
    nrow(annovar_input)
)
