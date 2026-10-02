# TrunKitten — reduced NMD-prediction annotation + scoring pipeline

End-to-end pipeline for annotating externally-called PTC / stop-gain variants with
the 8 features required by **TrunKitten**, the reduced 8-feature model derived from
**TrunCat** (TRUNcation-aware Classifier using Annotated Transcripts).

- **Input**: pre-called PTC variants (one transcript per variant via `txnames`).
- **Output**: 8-feature table + NMD escape predictions from TrunKitten.
- Feature conventions match the TrunCat training pipeline in
  `CobanAkdemirLab/NMDpredictionmodel` (see `DESIGN.md` for every ambiguity
  resolved + rationale).

> **Naming note**: the Python package is named `minicat` for historical reasons. All
> imports, module paths, and the `python -m minicat.cli` entry point continue to work
> unchanged. "TrunKitten" is the user-facing name for the reduced pipeline; `minicat`
> is its implementation package.

## The 8 features

`last.EJC`, `relativePTClocation`, `half_life_PC1`, `cdsseqs_AU_content`, `mut.exon`,
`cdsseqs_UC_content`, `phastcons_new3utr_first200_median`, `MedianExpression_log2`
(this is the model's column order).

## Install

```bash
pip install pandas numpy pyfaidx pyBigWig pyranges gffutils openpyxl catboost joblib pyyaml pytest
```

## Inputs

| File | Required columns / notes |
| --- | --- |
| `variants.tsv` | `variant_id, contig, position, refAllele, altAllele, gene, txnames` |
| `annotation.gtf(.gz)` | GENCODE v26 (hg38); `transcript_id`, `gene_id`, `exon_number` attributes |
| `genome.fa(.fai)` | hg38; must be indexed (`samtools faidx genome.fa`) |
| `phastcons.bw` | phastCons 100-way vertebrate, hg38; BigWig preferred, bedGraph also supported |
| `half_life_pc1.xlsx` | Agarwal & Kelley (2022) supplementary half-life table (column names set in the config) |
| `gtex_gene_median_tpm.gct.gz` | GTEx v8 "Median gene-level TPM by tissue" (`GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`), from <https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression>. The original `.gct` (two header lines) or a plain tab-separated table is accepted. |

phyloP is **not** needed any more.

Reference files live in `Model/TrunKitten/inputs/` (`../inputs/` from `pipeline/`; see `../inputs/README.md` for sources). All paths are set in `config/config.yaml`. Run the commands below from `pipeline/`.

## Run

```bash
# 1. Annotate
python -m minicat.cli --config config/config.yaml

# 2. Score
python predict.py \
    --annotated outputs/annotated.tsv \
    --model     ../model/trunkitten.pkl \
    --metadata  ../model/trunkitten_features.json \
    --out       outputs/predictions.tsv \
    --training-medians ../../TrunCat/predict/training_medians.json
```

Outputs:

- `outputs/annotated.tsv` — feature table, one row per input variant
- `outputs/qc_report.tsv` — exon counts, region lengths, boundary flags, missingness
  (`tx_not_in_gtf`, `half_life_missing`, `expression_missing`, ...)
- `outputs/run.log` — full annotation log
- `outputs/predictions.tsv` — TrunKitten escape probability and classification at the
  Youden threshold stored in `trunkitten_features.json`

## Feature definitions

See `DESIGN.md` §1. Every ambiguous feature (sequence windows, "after PTC"
conventions, DNA vs RNA alphabet, coding vs all exons) is addressed with a chosen
interpretation, rationale and code pointer. `MedianExpression_log2` follows
`Features/04_gene_level_features.R` (median across GTEx tissues of the per-tissue
median TPM, then `log2(x + 1)`).

## Important: feeding external cohorts to TrunKitten

The TrunCat training pipeline (notebook 02) imputes missing values as follows, and
`predict.py` applies the same rules:

- **Median-impute** `half_life_PC1` and `MedianExpression_log2` from the training
  medians. Pass `Model/TrunCat/predict/training_medians.json` via `--training-medians`
  (a flat `{column: median}` file shared with TrunCat's own predict notebooks).
- **No zero-fill** for any of the 8 features. A missing
  `phastcons_new3utr_first200_median` is left as NaN for CatBoost's native handling,
  exactly as in training and in TrunCat's own scoring.

If `--training-medians` is omitted, CatBoost's native NaN handling is used instead —
defensible, but it introduces a mild distribution shift.

## Validate against the training features

Check that the annotation pipeline reproduces the features in the TOPMed training
table (use the corrected merged table):

```bash
python scripts/validate_against_training.py \
    --merged ../../TrunCat/data/TOPMed_merged.csv \
    --config config/config.yaml \
    --n 300 --sample-seed 42 \
    --out outputs/validation_report.tsv
```

A random sample covers both strands, last-exon and upstream PTCs, and genes absent
from the half-life / GTEx tables. Add `--min-match 0.98` to make the script exit
non-zero if any feature matches in fewer than 98% of variants.

## Tests

```bash
pytest -v tests
```

`tests/test_feature_contract.py` checks that `minicat.REQUIRED_FEATURES` equals the
trained model's feature list and order, and unit-tests AU/UC content, the GTEx
expression table and the feature annotator on toy data.

## What this pipeline does NOT do

- **Call PTCs.** Inputs must already be confirmed stop-gained SNVs with a valid
  per-variant transcript choice in `txnames`.
- **Select transcripts.** TrunCat training-pipeline logic picks from multi-transcript
  `txnames` lists. TrunKitten expects exactly one transcript per row.
- **Validate the variant is a stop-gain.** We don't translate the CDS to verify
  ref→alt creates a PTC. That's upstream.
- **Run ensembles or calibration.** Single-model prediction at the training Youden
  threshold.

## On the relationship between TrunCat and TrunKitten

**TrunCat** is the full NMD-prediction model trained on 619 features. **TrunKitten**
is the reduced model trained on the 8 features that are in TrunCat's top 10 by mean
|SHAP| and in all five leave-one-fold-out top-10 rankings, intended for external-cohort
scoring where reproducing the full annotation pipeline is impractical. This repository
implements the TrunKitten annotation and scoring workflow.
