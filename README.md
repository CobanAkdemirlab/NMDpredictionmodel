[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/CobanAkdemirlab/NMDpredictionmodel/blob/main/Run_NMDpredictionmodel_pipeline.ipynb)
# NMDpredictionmodel

This repository contains the computational workflow used to annotate premature termination codon (PTC) variants, generate genomic and transcriptomic features, and develop machine-learning models for predicting nonsense-mediated mRNA decay (NMD) efficiency.

The workflow was developed using allele-specific expression (ASE) data from TOPMed and applied to independent population and disease datasets including gnomAD, ClinVar, and GREGoR.

---

## Overview

The analysis consists of four major stages:

1. **Dataset extraction and preprocessing**
2. **PTC and transcript annotation**
3. **Feature generation**
4. **Machine-learning model development and evaluation**

TOPMed was used as the model-development dataset because matched genomic and allele-specific expression information was available.

gnomAD, ClinVar, and GREGoR were processed using the same annotation and feature-generation framework for external application of the trained model.

---

## Pipeline Overview

```text
                         ┌── TOPMed ──────┐
                         ├── gnomAD ──────┤
Raw datasets ────────────├── ClinVar ─────┤
                         └── GREGoR ──────┘
                                  │
                                  ▼
                        Dataset extraction
                                  │
                                  ▼
                       Variant annotation
                                  │
                                  ▼
                 GENCODE v26 transcript structure
                                  │
                                  ▼
                         PTC annotation
                                  │
                                  ▼
                   Shared variant filtering
                                  │
                                  ▼
                       Feature generation
                                  │
                                  ▼
                      Feature matrix assembly
                                  │
                     ┌────────────┴────────────┐
                     ▼                         ▼
             TOPMed model training      External application
                                       gnomAD / ClinVar /
                                            GREGoR
```

For TOPMed, an additional ASE-processing step is used to derive the NMD outcome variable for model development.

---

## Repository Structure

```text
NMDpredictionmodel/
│
├── Dataset_extraction/
│   ├── 01a_TOPMed_preprocessing.sh
│   ├── 01b_TOPMed_extraction.R
│   ├── 02_gnomAD_extraction.R
│   ├── 03_ClinVar_extraction.R
│   ├── 04_GREGoR_extraction.R
│   └── README.md
│
├── Annotation/
│   ├── 01_variant_annotation.R
│   ├── 02_GENCODEv26_transcript_structure.R
│   ├── 03_PTC_annotation.R
│   ├── 04_shared_variant_filtering.R
│   ├── 05_TOPMed_ASE_simulation.R
│   └── README.md
│
├── Features/
│   ├── 01_GENCODEv26_sequence_features.R
│   ├── 02_PTBP1_binding_features.R
│   ├── 03_PTC_amino_acid_context.R
│   ├── 04_gene_level_features.R
│   ├── 05_PTC_geometry_derived_features.R
│   ├── 06a_prepare_variant_scores.R
│   ├── 06b_run_variant_scores.sh
│   ├── 06c_merge_variant_scores.R
│   ├── 07a_motif_region_extraction.py
│   ├── 07b_FIMO_motif_features.py
│   ├── 08a_make_optimal_codons_from_trna.py
│   ├── 08b_codon_optimality_features.py
│   ├── 09_conservation_score_features.py
│   ├── 10_EJC_occupancy_features.py
│   ├── 11_PTC_AUG_features.py
│   ├── 12_readthrough_features.py
│   ├── 99_build_feature_matrix.R
│   │
│   ├── helpers/
│   │   └── transcript_to_gene_mapping.py
│   │
│   ├── reference/
│   │   └── codon_optimality/
│   │       ├── hg38_UCSC_tRNA_table.tsv
│   │       └── optimal_codons.txt
│   │
│   └── README.md
│
├── Model/
│   ├── TrunCat/
│   └── TrunKitten/
│
├── Plotting/
│   └── [figure-generation scripts]
│
└── README.md
```

---

## 1. Dataset Extraction

Dataset-specific preprocessing scripts are located in `Dataset_extraction/`.

Four variant datasets are processed:

| Dataset | Role in analysis |
|---|---|
| **TOPMed Freeze 9b** | Model development and ASE-based NMD outcome |
| **gnomAD v4.1** | External population dataset |
| **ClinVar** | External clinically ascertained variant dataset |
| **GREGoR** | External rare-disease dataset |

Candidate protein-truncating variants were initially identified from the source datasets. Downstream analyses in this repository were restricted to **single-nucleotide PTC variants (SNVs)**.

Each dataset is converted to a standardized variant representation before entering the shared annotation and feature-generation workflow.

See `Dataset_extraction/README.md` for dataset-specific processing details.

---

## 2. Variant and PTC Annotation

Shared annotation scripts are located in `Annotation/`.

### Variant annotation

`01_variant_annotation.R`

Processes ANNOVAR transcript annotations and selects the transcript annotation used for downstream analyses.

### GENCODE transcript structure

`02_GENCODEv26_transcript_structure.R`

Constructs transcript-level reference information from **GENCODE v26**, including CDS structure, exon structure, CDS length, and exon count.

### PTC annotation

`03_PTC_annotation.R`

Maps each PTC to its coding transcript and derives core NMD-related positional annotations, including:

- coding position
- PTC-bearing exon
- distance from the PTC to the downstream exon junction
- distance from the PTC to the start and end of the CDS
- last-exon status
- penultimate-exon / last-50-nt status
- number of downstream exon junctions

### Shared variant filtering

`04_shared_variant_filtering.R`

Applies the common eligibility and filtering criteria used across datasets, including transcript structure and expression-related filters.

### TOPMed ASE outcome

`05_TOPMed_ASE_simulation.R`

TOPMed-specific processing is used to construct the ASE-derived outcome for model development.

The reference-allele ratio is calculated as:

```text
reference allele reads
────────────────────────────
reference + alternate reads
```

Repeatedly observed variants are handled using the carrier-sampling procedure implemented in the script.

---

## 3. Feature Generation

Feature-generation scripts are located in `Features/`.

Features span multiple biological levels.

### Transcript and sequence features

Derived from GENCODE v26 transcript structure and sequence, including:

- CDS length
- 5′ and 3′ UTR characteristics
- exon architecture
- nucleotide and dinucleotide composition
- transcript sequence context

### PTC geometry

PTC-specific positional features include:

- relative PTC position within the CDS
- PTC-bearing exon
- number of exons before and after the PTC
- distance from the PTC to exon boundaries
- relative position within the PTC-bearing exon

### Gene-level features

Gene-level annotations include:

- GTEx expression
- gnomAD constraint metrics
- mRNA half-life features

### Variant-level scores

Variant-level annotations include scores derived using ANNOVAR-compatible annotation databases, including:

- CADD
- REVEL
- gnomAD allele frequency

### RNA-binding protein motifs

Sequence regions surrounding the PTC, downstream exon junction, and untranslated regions are extracted and scanned using **FIMO (MEME Suite)**.

Motif matches are converted into variant-level RNA-binding protein feature matrices.

### Codon optimality

Optimal codons are derived from the **hg38 GtRNAdb-based tRNA annotation available through the UCSC Table Browser**.

`08a_make_optimal_codons_from_trna.py` reverse-complements each anticodon, counts tRNA gene copies per codon, and selects the most highly represented codon(s) for each amino acid.

The resulting `optimal_codons.txt` is used by `08b_codon_optimality_features.py` to calculate CDS-wide and PTC-local codon-optimality features.

### Conservation

Conservation features are calculated from the hg38:

- **phastCons100way**
- **phyloP100way**

BigWig tracks.

### EJC occupancy

Experimentally determined EJC occupancy information is derived from RIP-seq data associated with **NCBI GEO accession GSE41154**.

### Downstream AUG and Kozak context

Downstream AUG codons are identified in the reference transcript sequence in the original, +1, and +2 reading frames.

The workflow also evaluates sequence context surrounding downstream AUGs, including Kozak-context features.

### Readthrough features

PTC sequence context is evaluated using the readthrough scoring procedure implemented in `12_readthrough_features.py`.

---

## 4. Feature Matrix Assembly

`Features/99_build_feature_matrix.R` combines the individual feature tables using the standardized variant identifier.

The TOPMed model-development dataset produces the feature matrix used for machine-learning analysis:

```text
5,749 variants × 853 features
```

The same feature definitions are used when constructing matrices for gnomAD, ClinVar, and GREGoR.

Feature names and ordering are aligned to the TOPMed training feature set before external prediction.

---

## 5. Model Development and External Prediction

Model-development and prediction workflows are located in `Model/`.

### TrunCat

TrunCat is the primary classifier developed using the TOPMed ASE-derived NMD
outcome and the integrated genomic and transcriptomic feature matrix.

The TrunCat directory separates **model development using TOPMed** from
**application of the trained model to independent datasets**:

```text id="f85n7z"
Model/TrunCat/
│
├── data/
│   └── TOPMed model-development datasets
│
├── notebooks/
│   └── TrunCat data preparation, feature preprocessing,
│       model training, evaluation, and interpretation
│
├── model/
│   └── trained TrunCat model objects
│
├── results/
│   └── model-development results
│
└── predict/
    └── external prediction workflow and prediction files
```
---

## 6. Visualization and Exploratory Analyses

Scripts used for manuscript figures and exploratory analyses are maintained separately from the primary feature-generation workflow.

These analyses include:

- feature distributions
- NMD efficiency by biological feature
- continuous PTC-position relationships
- model performance
- SHAP feature importance
- external prediction summaries

Variables created only for visualization or exploratory categorization are not treated as primary feature annotations unless they were included in the model feature matrix.

---

## Software Requirements

### R

R ≥ 4.2

Major R packages used across the workflow include:

- `dplyr`
- `tidyr`
- `data.table`
- `GenomicRanges`
- `GenomicFeatures`
- `Biostrings`
- `rtracklayer`

### Python

Python ≥ 3.8

Major Python packages include:

- `pandas`
- `numpy`
- `pysam`
- `pyfaidx`
- `gffutils`
- `pyBigWig`
- `Biopython`

### External software

- ANNOVAR
- bcftools
- MEME Suite / FIMO

The motif-analysis workflow was performed using **MEME Suite 5.5.5**.

---

## Reference Genome and Transcript Annotation

Unless otherwise specified, analyses use:

- **Genome assembly:** GRCh38 / hg38
- **Transcript annotation:** GENCODE v26
- **ANNOVAR gene annotation:** `ensGene`

Additional feature-specific reference resources are documented in `Features/README.md` and within the corresponding scripts.

---

## Running the Workflow

The main workflow should be run in the following order:

```text
1. Dataset_extraction/
          ↓
2. Annotation/
          ↓
3. Features/
          ↓
   99_build_feature_matrix.R
          ↓
4. Model/
          ↓
5. Plotting/
```

Scripts within the feature-generation stage represent different feature families and are not necessarily sequential unless explicitly indicated.

For example:

```text
07a_motif_region_extraction.py
              ↓
07b_FIMO_motif_features.py
```

and:

```text
08a_make_optimal_codons_from_trna.py
              ↓
08b_codon_optimality_features.py
```

---

## Reproducibility

The repository has been organized so that each major analysis step has an explicit input and output and can be traced from the source variant datasets to the final model feature matrix.

Dataset-specific preprocessing is separated from shared annotation and feature-generation procedures.

Large source datasets that cannot be distributed through GitHub are documented with their source and version where possible. Smaller derived reference files required for reproducibility are included when redistribution permits.

Local file paths should be configured before running the workflow rather than being interpreted as part of the analysis specification.

---

## Outputs

The workflow produces:

- standardized PTC variant datasets
- transcript and PTC annotations
- individual biological feature tables
- the integrated model feature matrix
- trained NMD prediction models
- model performance results
- SHAP feature-importance results
- manuscript figures
