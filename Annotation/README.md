# Annotation

This directory contains the shared transcript and PTC annotation workflow
used for TOPMed, gnomAD, ClinVar, and GREGoR, together with the
TOPMed-specific ASE processing used for model development.

## Dataset usage

`01_variant_annotation.R` and `02_PTC_feature_annotation.R` are applied
to all datasets.

`03_TOPMed_ASE_simulation.R` is specific to TOPMed because TOPMed provides
the allele-specific expression outcome used for TrunCat model training
and internal evaluation.

gnomAD, ClinVar, and GREGoR proceed from the shared annotation steps
directly to feature generation and prediction.
