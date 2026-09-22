# TrunKitten — reduced NMD-prediction annotation + scoring pipeline

End-to-end pipeline for annotating externally-called PTC / stop-gain variants
with the 8 features required by **TrunKitten**, the reduced top-8 feature
model derived from **TrunCat** (TRUNcation-aware Classifier using Annotated
Transcripts). TrunKitten is the reduced TrunCat model using the top
informative transcript features.

- **Input**: pre-called PTC variants (one transcript per variant via `txnames`).
- **Output**: 8-feature table + NMD escape predictions from TrunKitten.
- Feature conventions match the TrunCat training pipeline in
  `CobanAkdemirLab/NMDpredictionmodel` exactly (see `DESIGN.md` for every
  ambiguity resolved + rationale).

> **Naming note**: the Python package is named `minicat` for historical
> reasons. All imports, module paths, and the `python -m minicat.cli`
> entry point continue to work unchanged. "TrunKitten" is the user-facing
> name for the reduced pipeline; `minicat` is its implementation package.

## Install

```bash
pip install pandas numpy pyfaidx pyBigWig pyranges gffutils openpyxl catboost joblib pyyaml pytest
```

## Inputs

| File                       | Required columns / notes                                                                     |
| -------------------------- | -------------------------------------------------------------------------------------------- |
| `variants.tsv`             | `variant_id, contig, position, refAllele, altAllele, gene, txnames`                          |
| `annotation.gtf(.gz)`      | Gencode-style; `transcript_id`, `gene_id`, `exon_number` attributes                          |
| `genome.fa(.fai)`          | Must be indexed (`samtools faidx genome.fa`)                                                 |
| `phastcons.bw`, `phylop.bw`| BigWig preferred; bedGraph also supported                                                    |
| `half_life_pc1.xlsx`       | -feature table (TrunKitten features), one row per input variant
- `outputs/qc_report.tsv` — exon counts, region lengths, boundary flags, missingness
- `outputs/run.log`       — full annotation log
- `outputs/predictions.tsv` — TrunKitten escape probability + classification at the training-Youden threshold

## Feature definitions

See `DESIGN.md` §1. Every ambiguous feature (sequence windows, "after PTC"
conventions, DNA vs RNA alphabet, coding vs all exons) is explicitly
addressed with a chosen interpretation + rationale + code pointer.

## Important: feeding external cohorts to TrunKitten

The TrunCat training pipeline applies specific imputation rules (Notebook 02):
- Median-impute: `half_life_PC1`, `MedianExpression_log2`, `CADD_phred`,
  `readthrough_score_hek293t` (of these, only `half_life_PC1` is in TrunKitten's
  top-8).
- Zero-fill: `phastcons_new3utr_first200_median`,
  `phylop_ptc_to_ejc_median` (when region is structurally absent).

Pass `Model/TrunCat/predict/training_medians.json` to `predict.py` via
`--training-medians`. It's a flat `{column: median}` file shared with
TrunCat's own predict notebooks; only the column relevant to TrunKitten
(`half_life_PC1`) is applied. Zero-fill for `phastcons_new3utr_first200_median`
and `phylop_ptc_to_ejc_median` is handled separately in `predict.py`
(`ZERO_FILL_COLUMNS`), since those aren't medians and aren't in that file.
If `--training-medians` is omitted, CatBoost's native NaN handling kicks
in — defensible, but introduces a mild distribution shift.

## What this pipeline does NOT do

- **Call PTCs.** Inputs must already be confirmed stop-gained SNVs with a
  valid per-variant transcript choice in `txnames`.
- **Select transcripts.** TrunCat training-pipeline logic picks from
  multi-transcript txnames lists. TrunKitten expects exactly one transcript
  per row — the assumption was called out in the spec.
- **Validate the variant is a stop-gain.** We don't translate the CDS to
  verify ref→alt creates a PTC. That's upstream.
- **Run ensembles or calibration.** Single-model prediction at the TrunCat
  training Youden threshold.

## On the relationship between TrunCat and TrunKitten

**TrunCat** is the full NMD-prediction model trained on all ~853 features.
**TrunKitten** is the reduced model trained on the top 8 features by mean
|SHAP| from TrunCat, intended for external-cohort scoring where reproducing
the full ~853-feature annotation pipeline is impractical. This repository
implements the TrunKitten annotation and scoring workflow.
