# TrunKitten — reduced NMD-prediction annotation + scoring pipeline

End-to-end pipeline for annotating externally-called PTC / stop-gain variants
with the 10 features required by **TrunKitten**, the reduced top-10 feature
model derived from **TrunCat** (TRUNcation-aware Classifier using Annotated
Transcripts). TrunKitten is the reduced TrunCat model using the top
informative transcript features.

- **Input**: pre-called PTC variants (one transcript per variant via `txnames`).
- **Output**: 10-feature table + NMD escape predictions from TrunKitten.
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
| `half_life_pc1.xlsx`       | Sheet with ENSG column and `half_life_PC1` column (configurable)                             |

Edit `config/config.yaml` to point to your files.

## Run

```bash
# 1. Annotate variants → 10-feature table + QC report (TrunKitten features)
python -m minicat.cli --config config/config.yaml

# 2. Score with the TrunKitten model
python predict.py \
    --annotated        outputs/annotated.tsv \
    --model            /path/to/models/trunkitten/trunkitten.pkl \
    --metadata         /path/to/results/trunkitten/top10_features.json \
    --training-medians /path/to/results/trunkitten/training_medians.json \
    --out              outputs/predictions.tsv
```

The `--model` and `--metadata` paths should point to the TrunKitten
artifacts produced by the reduced-model training notebook. The file names
shown above reflect the TrunKitten naming; if your saved artifacts still
use the historical `reduced_top10` naming, pass those paths — the scripts
don't care what the files are called, only what they contain.

## Validate

```bash
pytest -xvs tests/test_toy_transcripts.py
```

Hand-crafted toy transcripts exercise strand handling, exon-boundary PTCs,
last.EJC categorisation, transcript-position math, and region builders —
all offline, no reference files needed.

## Outputs

- `outputs/annotated.tsv` — 10-feature table (TrunKitten features), one row per input variant
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
  top-10).
- Zero-fill: `cdsseq_AUcontentlast200`, `phastcons_new3utr_first200_median`,
  `phylop_ptc_to_ejc_median` (when region is structurally absent).

Before scoring, write a `training_medians.json` (from Notebook 02's fitted
values) and pass it to `predict.py` via `--training-medians`. This reproduces
training-time preprocessing exactly. If omitted, CatBoost's native NaN
handling kicks in — defensible, but introduces a mild distribution shift.

Example `training_medians.json`:
```json
{
  "median_impute": { "half_life_PC1": 1.1303 },
  "zero_fill": [
    "cdsseq_AUcontentlast200",
    "phastcons_new3utr_first200_median",
    "phylop_ptc_to_ejc_median"
  ]
}
```

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

**TrunCat** is the full NMD-prediction model trained on all ~730 features.
**TrunKitten** is the reduced model trained on the top 10 features by mean
|SHAP| from TrunCat, intended for external-cohort scoring where reproducing
the full ~730-feature annotation pipeline is impractical. This repository
implements the TrunKitten annotation and scoring workflow.
