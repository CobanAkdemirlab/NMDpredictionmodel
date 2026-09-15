# Annotation

This directory contains the shared transcript and premature termination
codon (PTC) annotation workflow used for TOPMed, gnomAD, ClinVar, and
GREGoR, together with the TOPMed-specific allele-specific expression
(ASE) processing used for model development.

## Workflow overview

```text
Dataset-specific ANNOVAR output
            ↓
01_variant_annotation.R
            ↓
Canonical transcript selection
            ↓
02_GENCODEv26_transcript_structure.R
            ↓
Transcript structural information
(cds_length, cds_exons, exon_count, etc.)
            ↓
03_PTC_feature_annotation.R
            ↓
PTC positional and canonical NMD-rule annotation
            ↓
        ┌───────────────┐
        ↓               ↓
     TOPMed        gnomAD / ClinVar / GREGoR
        ↓               ↓
04_TOPMed_ASE_       Features/
simulation.R
        ↓
     Features/
