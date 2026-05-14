# NMDpredictionmodel

A comprehensive pipeline for annotating premature termination codon (PTC) variants and predicting nonsense-mediated decay (NMD) outcomes using genomic, transcript, and sequence-level features.

---

## Overview

This repository implements an end-to-end workflow that integrates large-scale genomic datasets (TOPMed, gnomAD, ClinVar, GREGoR) to:

* Annotate protein-truncating variants
* Extract biologically meaningful features
* Train predictive models of NMD efficiency
* Plot important features

---

## Pipeline Overview

```
Data → Cleaning → Annotation → Feature Engineering → Modeling → Visualization
```

---

## Directory Structure

```
NMDpredictionmodel/

├── 1_Data/
│   ├── TOPMed_extraction.R
│   ├── gnomAD_extraction.R
│   ├── ClinVar_extraction.R
│   └── GREGoR_extraction.R
│
├── 2_Features/
│   ├── canonical_transcript_annotation.R
│   ├── gencodev26_features.R
│   ├── variant_annotation.R
│   ├── transcript_features.R
│   ├── simulation.R
│   ├── ptc_features.R
│   ├── motif_regions_extraction.py
│   ├── motif_fimo_matrix.py
│   ├── ptc_amino_acid.R
│   ├── median_expression_mRNA_half_life.R
│   ├── transcript_to_gene_mapping.py
│   ├── gnomAD_constraint.R
│   ├── ejc_occupancy.py
│   ├── ptc_aug_analyzer.py
│   ├── conservation_scores.py
│   ├── readthrough_scoring.py
│   ├── ejc_analysis.py
│   └── annotated_output/
│
├── 3_Model/
│   └── (ENTER HERE)
│
├── 4_Plotting/
│   ├── plot_predictions.py
│   ├── feature_distributions.py
│   └── figures/
│
└── README.md
```

---

## Key Components

### Data Processing

Extraction and cleaning of variant datasets from:

* TOPMed
* gnomAD
* ClinVar
* GREGoR

---

### Feature Engineering

Includes:

* Transcript structure (exon count, isoform count, CDS length)
* PTC positional features (distance to exon junctions, stop codon)
* Sequence composition (UTR, CDS, dinucleotide content)
* RNA-binding protein motifs (FIMO)
* Gene-level features (expression, constraint metrics)
* Conservation scores (PhastCons, PhyloP)

---

### Modeling

* Machine learning models trained to predict NMD escape
* Uses engineered features across multiple biological layers


---

### Visualization

* Feature distributions
* Model predictions
* SHAP plots

---

## Getting Started

### Requirements

* R (≥ 4.2)
* Python (≥ 3.8)

### External tools

* ANNOVAR
* bcftools
* MEME Suite (FIMO)

---

## Workflow

Run modules in order:

1. `1_Data/` → dataset extraction
2. `2_Features/` → annotation & feature generation
3. `3_Model/` → model training
5. `4_Plotting/` → visualization

---

## Output

* Annotated variant tables
* Feature matrices
* Trained models
* SHAP values
* Publication-ready figures

---

## Notes

* Paths may need adjustment for local environments
* Large datasets are not included
* Scripts are modular and can be run independently

