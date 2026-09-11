# Dataset Extraction

This directory contains the dataset-specific preprocessing workflows used to prepare variants for the TrunCat annotation and feature-generation pipeline.

## Overview

Four genomic datasets were processed in this study. **TOPMed Freeze 9b was used for TrunCat model development and internal evaluation**, including construction of the training feature matrix and model training/testing.

**gnomAD, ClinVar, and GREGoR were not used for model training.** These datasets were processed independently and used for application of the trained TrunCat model and generation of NMD predictions.

| Dataset | Version / release | Genome build | Variant type | Role |
|---|---|---|---|---|
| TOPMed | Freeze 9b | GRCh38 | SNVs | Model development, training, and internal testing |
| gnomAD | v4.1 | GRCh38 | SNVs | External prediction |
| ClinVar | 2026-02-01 | GRCh38 | SNVs | External prediction |
| GREGoR | [ADD DATA RELEASE] | GRCh38 | SNVs | External prediction |

> **Reproducing the TrunCat model-development dataset:** The TOPMed workflow (`01a` → `01b`) is the dataset-extraction path required to reconstruct the model-development feature matrix. gnomAD, ClinVar, and GREGoR are processed independently for external prediction analyses.

---

## Directory structure

```text
Dataset_extraction/
├── README.md
├── 01a_TOPMed_preprocessing.sh
├── 01b_TOPMed_extraction.R
├── 02_gnomAD_extraction.R
├── 03_ClinVar_extraction.R
└── 04_GREGoR_extraction.R
```

The general workflow is:

```text
Source variant data
        ↓
Dataset-specific preprocessing
        ↓
PTC/PTV identification
        ↓
SNV restriction
        ↓
Standardized GRCh38 variant representation
        ↓
ANNOVAR ensGene annotation
        ↓
../Annotation/
        ↓
../Features/
```

The datasets are processed independently. Numbering is provided to make the repository easier to navigate and does not indicate that gnomAD, ClinVar, or GREGoR depend on one another.

---

# Variant inclusion

The downstream analyses described in this repository were restricted to **single-nucleotide variants (SNVs)**.

Candidate protein-truncating or premature termination codon (PTC) variants were identified within each source dataset, followed by restriction to SNVs before downstream ANNOVAR annotation and feature generation.

This provides a consistent SNV-based variant representation across TOPMed, gnomAD, ClinVar, and GREGoR.

Unless otherwise specified, downstream PTC analyses in this workflow therefore refer to the SNV subset retained for analysis.

---

# 1. TOPMed

## Role in TrunCat

TOPMed Freeze 9b was used as the **model-development dataset** for TrunCat.

The TOPMed workflow generates the allele-specific expression (ASE) and genomic variant information that ultimately contributes to construction of the model feature matrix used for model training and internal testing.

TOPMed individual-level genomic and transcriptomic data are controlled-access and therefore cannot be distributed through this repository.

## Workflow

TOPMed preprocessing consists of two sequential scripts:

```text
01a_TOPMed_preprocessing.sh
            ↓
merged_TOR*.vcf
            ↓
01b_TOPMed_extraction.R
            ↓
../Annotation/
```

---

## Step 1 — Candidate PTV extraction, ASE quantification, and genotype matching

Run:

```bash
bash 01a_TOPMed_preprocessing.sh
```

### Required input

The preprocessing workflow requires:

- TOPMed Freeze 9b annotated variant data
- TOPMed Freeze 9b chromosome-specific genotype BCF files
- matched RNA-seq BAM files
- sample metadata linking RNA-seq TOR IDs to TOPMed genomic NWD IDs
- GRCh38 reference FASTA

### Processing

`01a_TOPMed_preprocessing.sh` performs the following steps:

1. extracts candidate variants annotated as `stop_gained` or `frameshift_variant`;
2. removes duplicate variant records;
3. indexes the resulting candidate PTV VCF;
4. runs GATK `ASEReadCounter` on matched RNA-seq BAM files;
5. identifies genomic positions with allele-specific RNA-seq measurements;
6. retrieves the corresponding donor genotype from the TOPMed Freeze 9b genotype BCFs;
7. matches RNA-seq ASE measurements to the corresponding donor genotype.

Candidate stop-gained and frameshift variants are initially considered during ASE extraction. The downstream analysis is subsequently restricted to SNVs in `01b_TOPMed_extraction.R`.

### Primary output

```text
TOPMed_PTVs.dedup.vcf.gz
ASE_results/*_ASE_counts.table
ASE_genotype/merged_TOR*.vcf
```

Each `merged_TOR*.vcf` file contains allele-specific RNA-seq measurements together with the corresponding donor genotype.

---

## Step 2 — Heterozygous SNV extraction and ANNOVAR preparation

Run:

```bash
Rscript 01b_TOPMed_extraction.R
```

### Input

```text
ASE_genotype/merged_TOR*.vcf
```

### Processing

`01b_TOPMed_extraction.R`:

1. reads the per-individual ASE/genotype files generated in Step 1;
2. retains heterozygous variant observations;
3. combines qualifying observations across individuals;
4. restricts the downstream analysis to SNVs;
5. constructs standardized genomic variant identifiers;
6. prepares variants for ANNOVAR;
7. performs transcript-level annotation using the hg38 `ensGene` database.

### Output

Representative outputs include:

```text
TOPMed_heterozygous_PTC_SNVs.tsv
TOPMed_PTC_SNV.vcf
TOPMed_PTC_SNV.avinput
TOPMed_PTC_SNV_ensGene.variant_function
TOPMed_PTC_SNV_ensGene.exonic_variant_function
```

These files are passed to the downstream annotation and feature-generation workflow.

---

# 2. gnomAD

## Role in TrunCat

gnomAD was **not used for TrunCat model training**. gnomAD variants were processed through the annotation and feature-generation pipeline for application of the trained TrunCat model.

## Data source

- Dataset: Genome Aggregation Database (gnomAD)
- Version: **v4.1**
- Genome build: **GRCh38**

Initial download and chromosome-level preprocessing were performed using the gnomAD extraction workflow described here:

[ADD EXACT LINK TO `gnomAD_downloaddata.R`]

The resulting chromosome-level RDS objects are used as input to:

```text
02_gnomAD_extraction.R
```

## Processing

`02_gnomAD_extraction.R`:

1. reads and merges chromosome-level gnomAD variant objects;
2. processes variants using `aenmd`;
3. removes variants containing undefined alternative alleles;
4. identifies variants predicted to introduce a premature termination codon;
5. restricts the dataset to SNVs;
6. creates standardized genomic variant identifiers;
7. generates ANNOVAR input;
8. performs transcript-level annotation using the hg38 `ensGene` database.

## Output

Representative outputs include:

```text
gnomAD_v4.1_all_variants.rds
gnomAD_v4.1_PTC_SNV.rds
gnomAD_v4.1_PTC_SNV.vcf
gnomAD_v4.1_PTC_SNV.avinput
gnomAD_v4.1_PTC_SNV_ensGene_*
gnomAD_v4.1_fr.var.can_snv.RData
```

These files are subsequently processed by the downstream annotation and feature-generation workflow.

---

# 3. ClinVar

## Role in TrunCat

ClinVar was **not used for TrunCat model training**. ClinVar variants were processed independently for application of the trained TrunCat model to clinically interpreted variants.

## Data source

- Dataset: NCBI ClinVar
- Release: **February 1, 2026**
- Input VCF: `clinvar_20260201.vcf.gz`
- Genome build: **GRCh38**

## Processing

Run:

```bash
Rscript 03_ClinVar_extraction.R
```

The script:

1. reads the ClinVar VCF;
2. processes variants using `aenmd`;
3. removes variants containing undefined alternative alleles;
4. annotates variants for predicted premature termination codons;
5. retains PTC-producing SNVs;
6. creates standardized genomic variant identifiers;
7. generates ANNOVAR input;
8. performs transcript-level annotation using the hg38 `ensGene` database.

## Output

Representative outputs include:

```text
clinvar_20260201_PTC_SNV.rds
clinvar_20260201_PTC_SNV.vcf
clinvar_20260201_PTC_SNV.avinput
clinvar_20260201_PTC_SNV_ensGene.variant_function
clinvar_20260201_PTC_SNV_ensGene.exonic_variant_function
```

These files are subsequently processed by the downstream annotation and feature-generation workflow.

---

# 4. GREGoR

## Role in TrunCat

GREGoR was **not used for TrunCat model training**. GREGoR variants were processed independently for application of the trained TrunCat model.

Individual-level GREGoR data are controlled-access and therefore cannot be distributed through this repository.

## Input

The extraction workflow begins with the deduplicated GREGoR variant file:

```text
filtered_oc_base_matches.unique.vcf
```

The input variants use GRCh38 genomic coordinates.

## Processing

Run:

```bash
Rscript 04_GREGoR_extraction.R
```

The script:

1. reads the GREGoR variant file;
2. standardizes chromosome naming for compatibility with `aenmd`;
3. processes variants using `aenmd`;
4. removes variants containing undefined alternative alleles;
5. identifies variants predicted to introduce a premature termination codon;
6. restricts the dataset to SNVs;
7. creates standardized genomic variant identifiers;
8. generates ANNOVAR input;
9. performs transcript-level annotation using the hg38 `ensGene` database.

## Output

Representative outputs include:

```text
GREGoR_PTC_SNV.rds
GREGoR_PTC_SNV.vcf
GREGoR_PTC_SNV.avinput
GREGoR_PTC_SNV_ensGene.variant_function
GREGoR_PTC_SNV_ensGene.exonic_variant_function
```

These files are subsequently processed by the downstream annotation and feature-generation workflow.

---

# Reference annotation

Unless otherwise specified, genomic coordinates and transcript annotations in this workflow use:

| Resource | Version |
|---|---|
| Genome assembly | GRCh38 / hg38 |
| GENCODE | v26 |
| ANNOVAR gene annotation | `ensGene` |

The same reference framework is maintained across datasets to support consistent downstream feature generation.

---

# Software requirements

The dataset-extraction workflow requires a combination of command-line and R software.

## Command-line software

The TOPMed preprocessing workflow requires:

```text
bcftools
tabix
GATK
```

ANNOVAR is required for downstream transcript-level annotation.

## R

Major R/Bioconductor dependencies include:

```text
aenmd
Biostrings
GenomicRanges
GenomeInfoDb
S4Vectors
VariantAnnotation
data.table
dplyr
stringr
```

Additional package requirements are documented within individual scripts.

Package versions used for the manuscript analysis should also be recorded in the repository-level software environment documentation.

---

# Configuring file paths

The original analyses were performed in a high-performance computing environment. User-specific absolute paths should not be required to reproduce the public workflow.

Where applicable, input directories, output directories, reference files, and software locations are defined in configuration sections rather than throughout the analysis code.

Users should modify these configuration values for their local computing environment before running the scripts.

For example:

```r
CONFIG <- list(
    input_dir   = "/path/to/input",
    output_dir  = "/path/to/output",
    annovar_dir = "/path/to/annovar"
)
```

---

# Reproducibility and expected workflow

The complete TrunCat model-development workflow proceeds from TOPMed through annotation and feature generation:

```text
TOPMed Freeze 9b
       ↓
Dataset_extraction/
       ↓
ASE + heterozygous PTC SNVs
       ↓
Annotation/
       ↓
Transcript-specific PTC annotations
       ↓
Features/
       ↓
Final model feature matrix
       ↓
Model/
       ↓
TrunCat
```

The model-development feature matrix reported in the study contains:

```text
5,749 variants × 853 features
```

gnomAD, ClinVar, and GREGoR do **not** contribute observations to model training. After TrunCat model development, these datasets are independently processed through the compatible annotation and feature-generation workflow for prediction.

See the repository-level README and the README files within `Annotation/`, `Features/`, and `Model/` for subsequent steps.
