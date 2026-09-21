# ==============================================================================
# Exploratory feature transformations
#
# Purpose:
#   Generate log2-transformed and categorized versions of selected transcript
#   length features for descriptive analyses, single-feature association tests,
#   and plotting.
#
# These variables were not used as primary transcript annotations and should
# not be confused with the continuous features used in the feature-generation
# workflow.
#
# Input variables:
#   - cds_length
#   - threeUTR_length
#   - fiveutr_length
#   - log2newUTR
#
# Derived exploratory variables:
#   - log2_CDS
#   - cds_length.cut
#   - log2_3utr
#   - threeUTR_length.cut
#   - log2_5utr
#   - fiveUTR_length.cut
#   - log2newUTR.cut
# ==============================================================================



