# External Validation Using LMNA Saturation Genome Editing (SGE) Data

This directory contains the external-validation workflow used to evaluate TrunKitten predictions against experimentally measured nonsense-mediated mRNA decay (NMD) efficiency from the *LMNA* saturation genome editing (SGE) dataset reported by Cortazar et al.

The analysis applies the pretrained **eight-feature TrunKitten model** to independently generated LMNA variants without retraining the model. Predictions are compared with normalized experimental NMD measurements to assess concordance between predicted and observed NMD efficiency.

## Workflow Overview

```text
LMNA SGE experimental data
           ↓
01_LMNA_prepare_TrunKitten_features.R
           ↓
SNV-compatible stop-gain variants
           ↓
Conservation annotation
           ↓
02_LMNA_add_conservation_and_predict.R
           ↓
Eight-feature TrunKitten input matrix
           ↓
Model/TrunKitten/pipeline/predict.py
           ↓
TrunKitten NMD predictions
           ↓
03_LMNA_SGE_TrunKitten_validation.R
           ↓
Experimental validation and visualization
```

## 1. Prepare LMNA Variants and Features

The first script processes the LMNA SGE dataset, identifies single-nucleotide stop-gain variants, and reconstructs transcript-level and sequence-based features using GENCODE v26 and GRCh38.

Run from the repository root:

```bash
Rscript External_validation/01_LMNA_prepare_TrunKitten_features.R \
  /path/to/NMD_efficiency_values_PTC.txt \
  /path/to/gencode.v26.primary_assembly.annotation.gtf \
  /path/to/GRCh38.primary_assembly.genome.fa \
  Model/TrunKitten/inputs/half_life_pc1.xlsx \
  Features/reference/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.txt \
  output
```

The script generates:

- `output/LMNA_Cortazar_TrunKitten_variants.csv`
- `output/lmna_snv_before_conservation.rds`

## 2. Add Conservation Features

The LMNA variants are annotated using the existing conservation feature-generation workflow.

The required conservation input is:

`output/variants_with_conservation_medians.csv`

This file must contain the variant identifier (`variant_id`) and the model-required conservation feature (`phastcons_new3utr_first200_median`).

Generate the final eight-feature input matrix:

```bash
Rscript External_validation/02_LMNA_add_conservation_and_predict.R \
  output/lmna_snv_before_conservation.rds \
  output/variants_with_conservation_medians.csv \
  output/LMNA_TrunKitten_8features_with_id.tsv
```

The resulting matrix contains eight model features plus the variant identifier.

## 3. Generate TrunKitten Predictions

Apply the pretrained TrunKitten model:

```bash
python Model/TrunKitten/pipeline/predict.py \
  --annotated output/LMNA_TrunKitten_8features_with_id.tsv \
  --model Model/TrunKitten/model/trunkitten.pkl \
  --metadata Model/TrunKitten/model/trunkitten_features.json \
  --out output/LMNA_TrunKitten_predictions.tsv
```

### Model Inputs

The eight features used by the current TrunKitten model are:

| Feature | Description |
|---|---|
| `last.EJC` | PTC position relative to the last exon junction |
| `relativePTClocation` | Relative PTC position within the coding sequence |
| `half_life_PC1` | mRNA half-life principal component |
| `cdsseqs_AU_content` | Coding-sequence AU content |
| `mut.exon` | PTC-bearing exon |
| `cdsseqs_UC_content` | Coding-sequence UC content |
| `phastcons_new3utr_first200_median` | Conservation score in the first 200 nt of the new 3′ UTR |
| `MedianExpression_log2` | Log2-transformed median GTEx expression |

The model metadata specifies a Youden classification threshold of **0.4313564920894303**. The pretrained model should be checked for compatibility with the supplied metadata before prediction.

### Prediction Output

`output/LMNA_TrunKitten_predictions.tsv`

This file contains variant-level TrunKitten predictions, including predicted NMD escape probabilities.

## 4. External Validation Against Experimental NMD

The final script compares TrunKitten predictions with the experimental LMNA SGE measurements.

Experimental PTC measurements are normalized using regional SNV measurements, following the normalization procedure implemented in the validation script.

Run:

```bash
Rscript External_validation/03_LMNA_SGE_TrunKitten_validation.R \
  /path/to/NMD_efficiency_values_PTC.txt \
  /path/to/NMD_efficiency_values_SNV.txt \
  output/lmna_snv_before_conservation.rds \
  output/LMNA_TrunKitten_predictions.tsv \
  output/validation
```

### Validation Analysis

The primary validation analysis calculates the **Spearman correlation** between:

- TrunKitten-predicted NMD probability (`1 - escape_prob`)
- Normalized experimental NMD efficiency

Variants that do not meet the experimental quality-control criterion (`std_val ≤ 0.2`) are excluded from the primary validation analysis.

The validation script generates the final comparison dataset, correlation statistics, and visualization.

### Expected Outputs

Results are written to `output/validation/`, including the final validation dataset:

`LMNA_SGE_validation_final.csv`
## Supplementary Figure S13 — LMNA SGE Validation

Supplementary Figure S13 evaluates the concordance between
eight-feature TrunKitten predictions and experimentally measured
NMD efficiency from the LMNA saturation genome editing dataset.

- **Panel A:** Spearman correlation between predicted and experimental
  NMD efficiency.
- **Panel B:** Comparison of experimental NMD efficiency between
  variants predicted to trigger or escape NMD.

Figure-generation script:

[04_LMNA_SGE_Supplementary_Figure_S13.R](04_LMNA_SGE_Supplementary_Figure_S13.R)

The figure includes 84 SNV-compatible stop-gain variants passing
experimental quality control.

## Reproducibility Notes

The LMNA SGE dataset is used for independent validation and is not included in TrunKitten model training.

The previous ten-feature workflow has been replaced by the eight-feature model configuration. Results from the earlier ten-feature analysis should not be interpreted as validation results for the updated model.

The historical workflow identified 86 SNV-compatible LMNA stop-gain variants and 84 experimental observations passing quality control. These values should be checked when reproducing the analysis.

All reported validation statistics should be generated from the current eight-feature model predictions.
