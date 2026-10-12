# NMDpredictionmodel

This repository contains the computational workflow for annotating premature termination codon (PTC) variants, generating genomic and transcriptomic features, and developing machine-learning models to predict nonsense-mediated mRNA decay (NMD) efficiency.

TrunCat was developed using allele-specific expression (ASE) data from **TOPMed Freeze 9b**. The trained model was subsequently applied to **gnomAD v4.1, ClinVar (February 1, 2026), and GREGoR**. These external datasets were not used for model training.

## Getting started with Google Colab

We provide a step-by-step ClinVar example notebook to demonstrate variant annotation, feature generation, and application of the trained TrunCat model.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/CobanAkdemirlab/NMDpredictionmodel/blob/main/TrunCat_ClinVar_End_to_End_GitHub.ipynb)

**[Open the TrunCat ClinVar tutorial](TrunCat_ClinVar_End_to_End_GitHub.ipynb)** · **[Example VCF](examples/)**

The tutorial is organized around five stages: environment setup; variant and PTC annotation; feature generation; model prediction; and results, including ClinVar clinical significance and visualization.

**Reproducibility status:** The notebook has not yet been verified to execute uninterrupted in a fresh Colab runtime. Some steps may require additional reference resources or configuration. The repository's scripts and directory READMEs document the individual processing stages.

## Pipeline overview

```text
TOPMed Freeze 9b ────────┐
gnomAD v4.1 ─────────────┤
ClinVar (2026-02-01) ────┼──> Dataset-specific extraction
GREGoR ──────────────────┘               |
                                        v
                              Variant / PTC annotation
                              (ANNOVAR, GENCODE v26)
                                        |
                                        v
                                Feature generation
                                        |
                                        v
                               Feature matrix assembly
                                        |
                           +------------+------------+
                           |                         |
                           v                         v
                     TOPMed training          External prediction
                     and evaluation        (gnomAD / ClinVar / GREGoR)
```

TOPMed additionally uses matched genomic and RNA-seq information to derive the ASE-based outcome for model development. The downstream variant analyses described here focus on PTC-producing single-nucleotide variants (SNVs).

## Workflow documentation

| Directory | What it contains |
|---|---|
| [`Dataset_extraction/`](Dataset_extraction/README.md) | Dataset-specific preparation of TOPMed, gnomAD, ClinVar, and GREGoR variants; ANNOVAR input preparation |
| [`Annotation/`](Annotation/README.md) | Shared transcript structure, PTC geometry and NMD-rule annotations; TOPMed ASE processing |
| [`Features/`](Features/README.md) | Sequence, transcript, gene-level, motif, conservation, codon, EJC, AUG, and readthrough features; feature assembly |
| [`Model/TrunCat/`](Model/TrunCat/) | Model development, trained CatBoost classifier, and prediction utilities |
| [`Model/TrunKitten/`](Model/TrunKitten/) | Additional model analyses |
| [`Plotting/`](Plotting/) | Manuscript figures and exploratory analyses |
| [`examples/`](examples/) | Example input files for testing the workflow |

For the required inputs, individual scripts, reference files, and processing details, follow the README in each directory.

## TrunCat prediction

The trained TrunCat classifier outputs a **predicted probability of NMD escape** (class 1). The documented decision threshold is **0.4747**:

- **NMD:** escape probability below 0.4747.
- **NMD escape:** escape probability at or above 0.4747.

The trained model and prediction utilities are under [`Model/TrunCat/`](Model/TrunCat/). ClinVar pathogenicity labels are separate clinical annotations, **not** experimentally observed NMD outcomes or model predictors.

## Reference framework and software

- **Genome:** GRCh38 / hg38
- **Transcript annotation:** GENCODE v26
- **ANNOVAR annotation:** `ensGene`
- **Main languages:** R and Python
- **External tools used in the workflow:** ANNOVAR, bcftools, MEME Suite / FIMO (motif analyses used MEME Suite 5.5.5)

Additional R/Bioconductor and Python dependencies, reference resources, and script-specific setup are described in the directory READMEs and source code. Controlled-access TOPMed and GREGoR data cannot be redistributed in this repository.

## Outputs

Depending on the workflow, outputs include standardized PTC variants, transcript/PTC annotations, biological feature tables, integrated model inputs, NMD predictions, model evaluation summaries, and figures. The TOPMed model-development analysis used **5,749 variants and 619 selected features**.
