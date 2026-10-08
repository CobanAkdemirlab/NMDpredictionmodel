# TrunCat

**TRUN**cation-aware **C**lassifier using **A**nnotated **T**ranscripts.

TrunCat is the full NMD-escape prediction pipeline: a CatBoost classifier
trained on 5,749 premature termination codon (PTC) variants from TOPMed
(3,214 genes; 45.5% NMD-escape), using 619 transcript-level, genomic, and
sequence-based features to predict whether stop-gain variants escape
nonsense-mediated decay (NMD).

**Out-of-fold ROC-AUC = 0.777** (5-fold gene-grouped CV; fold mean 0.779 ± 0.015) |
**PR-AUC = 0.733** | **Youden-optimal threshold = 0.47**

## Layout

```text
Model/TrunCat/
├── config/
│   ├── config.yaml                  paths, hyperparameters, feature lists
│   └── feature_column_order.txt     fixed column order for the merged table
├── data/                            TOPMed-derived input and intermediate tables
├── notebooks/
│   ├── 01_data_loading_and_merging.ipynb        pipeline stage 1
│   ├── 02_feature_cleaning_and_selection.ipynb  pipeline stage 2
│   ├── 03_model_training.ipynb                  pipeline stage 3
│   └── ...                          analysis notebooks (see "Analysis notebooks")
├── scripts/                         standalone Python equivalents of 01-03,
│                                    check_parity.py, requirements.txt, setup.py
├── model/                           trained TrunCat artifacts
│   ├── truncat.cbm
│   └── truncat.pkl
├── predict/                         scoring new cohorts (predict.py,
│                                    predict_cohorts.ipynb, training_medians.json)
├── track/                           UCSC track hub and bigBed of predictions
└── results/                         CV predictions, SHAP, ablation, LOFO,
                                     calibration and tier outputs
```

## Pipeline

Notebooks and scripts share `config/config.yaml` and produce the same outputs;
`scripts/check_parity.py` verifies this for stages 01-03. Use whichever
interface fits your workflow.

1. **01** - merge TOPMed variants with six annotation sources into `data/TOPMed_merged.csv`, with columns put in the order given by `config/feature_column_order.txt`
2. **02** - feature cleaning, leakage removal, correlation filtering -> `data/TOPMed_cleaned.csv`
3. **03** - 5-fold gene-grouped stratified cross-validation, SHAP analysis, final model -> `model/truncat.{cbm,pkl}` and `results/`

See [`notebooks/README.md`](notebooks/README.md) for a per-notebook walkthrough
and [`scripts/README.md`](scripts/README.md) for the script interface.

### Reproducibility notes

- **Fixed iterations, no early stopping.** Every fold model and the final model use the same fixed number of boosting iterations (250), so out-of-fold performance is not tuned on the evaluation fold.
- **Gene-grouped folds.** `StratifiedGroupKFold` (seed 42) keeps all variants of a gene in one fold.
- **scikit-learn 1.3-1.7.** scikit-learn 1.8 changes the `StratifiedGroupKFold` splits, so the committed folds are not reproduced with it. Tested with Python 3.12, scikit-learn 1.5.1, pandas 2.3.3, catboost 1.2.8, numpy 1.26.4, shap 0.48.0.
- **Column order matters.** CatBoost results shift by about 0.002 AUC if the column order changes; keep `config/feature_column_order.txt`.

## Features

Mean |SHAP| rankings are in `results/feature_importances_cv_averaged.csv`.
Eight features are in the global top 10 *and* in every leave-one-fold-out top 10
(`results/lofo_stability/`): `last.EJC`, `relativePTClocation`, `half_life_PC1`,
`cdsseqs_AU_content`, `mut.exon`, `cdsseqs_UC_content`,
`phastcons_new3utr_first200_median`, `MedianExpression_log2`. These define
TrunKitten.

## Analysis notebooks

Run from `notebooks/`. Outputs go to `results/`.

| Notebook | Purpose |
|---|---|
| `add_ids_to_cv_predictions.ipynb` | attach variant and gene IDs to the OOF predictions |
| `lofo_stability_analysis.ipynb` | leave-one-fold-out stability of the SHAP top 10 |
| `ablation_analysis.ipynb` | feature-removal curve, redundancy, 50-nt rule ablation, TrunKitten variants |
| `trunkitten_bootstrap_ci.ipynb` | gene-clustered paired bootstrap, TrunKitten vs TrunCat |
| `calibration_analysis.ipynb` | calibration metrics and confidence-tier check on OOF predictions; writes `results/calibration/tier_thresholds.json` |
| `tier_analysis_clinvar_combined.ipynb` | tier composition of ClinVar, gnomAD and GREGoR predictions (run after the calibration notebook) |
| `build_ucsc_track.ipynb` | UCSC track hub from cohort predictions |

## Scoring new variants

`predict/predict_cohorts.ipynb` scores the ClinVar, gnomAD and GREGoR tables;
`predict/predict.py` is the scoring function. Input tables must contain the
annotation columns the model was trained on (listed in the `feature` column of
`data/final_feature_list.csv`), generated with the scripts in `Features/`.

## UCSC Genome Browser track

TrunCat predictions for the scored gnomAD, ClinVar and GREGoR variants are
available as a UCSC track hub (hg38):
[load in the Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTracks?hubUrl=https://raw.githubusercontent.com/CobanAkdemirLab/NMDpredictionmodel/main/Model/TrunCat/track/hub.txt&db=hg38).
Files and rebuild instructions are in [`track/`](track/README.md).

## Related

For a lightweight version using only 8 features - suitable for scoring
external cohorts without reproducing the full annotation pipeline - see
[`../TrunKitten/`](../TrunKitten/) (OOF ROC-AUC = 0.775).
