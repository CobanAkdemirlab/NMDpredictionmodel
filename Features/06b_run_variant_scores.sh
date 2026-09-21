#!/usr/bin/env bash

# ==============================================================================
# ANNOVAR variant-level score annotation
#
# Purpose:
#   Annotate PTC SNVs with variant-level pathogenicity and population
#   frequency information.
#
# Genome build:
#   GRCh38 / hg38
#
# ANNOVAR databases used in the manuscript analysis:
#
#   dbnsfp42a
#       dbNSFP version 4.2a
#       ANNOVAR hg38 release: 2021-07-10
#
#   gnomad_exome
#       gnomAD exome collection version 2.0.1
#
# IMPORTANT:
#   These database versions are intentionally fixed to reproduce the
#   feature values used in the manuscript. Do not substitute newer
#   ANNOVAR databases when reproducing the published feature matrix.
#
# Input:
#   variants_for_external_annotation.avinput
#
# Output:
#   variants_external_annotation.hg38_multianno.txt
# ==============================================================================

set -euo pipefail


# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

ANNOVAR_DIR="/path/to/annovar"
HUMANDB="${ANNOVAR_DIR}/humandb"

INPUT_FILE="variants_for_external_annotation.avinput"
OUTPUT_PREFIX="variants_external_annotation"


# ------------------------------------------------------------------------------
# Run ANNOVAR
# ------------------------------------------------------------------------------

perl "${ANNOVAR_DIR}/table_annovar.pl" \
    "${INPUT_FILE}" \
    "${HUMANDB}" \
    -buildver hg38 \
    -out "${OUTPUT_PREFIX}" \
    -remove \
    -protocol dbnsfp42a,gnomad_exome \
    -operation f,f \
    -nastring . \
    -otherinfo

echo "Variant-level ANNOVAR annotation complete."
