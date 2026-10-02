# TrunCat — Standalone Scripts

Python script equivalents of the three notebooks. Same inputs, same outputs,
same config. Use when you want command-line execution instead of Jupyter, or
when another notebook (for example `Run_NMDpredictionmodel_pipeline.ipynb`)
calls the pipeline.

## Requirements

- Python 3.12+
- Dependencies in `requirements.txt`

```bash
python -m venv venv
source venv/bin/activate       # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

The committed results were produced with Python 3.12.2, scikit-learn 1.5.1, pandas 2.3.3,
CatBoost 1.2.8, numpy 1.26.4 and shap 0.48.0.

**scikit-learn must be below 1.8.** Version 1.8 builds different gene-grouped
cross-validation folds, so the reported CV numbers cannot be reproduced with it.
The reported results were checked with scikit-learn 1.3 to 1.7. pandas 2.x and
3.x both give identical cleaned data. CatBoost results vary slightly with the
CatBoost version and thread count (out-of-fold AUC within about 0.001 in our
checks), so expect close agreement, not bit-identical models.

## Verify setup

```bash
python setup.py
```

Checks config validity, input file presence, output directory creation, and
dependency installation (including the scikit-learn version above).

## Usage

Run in order from `Model/TrunCat/scripts/`:

```bash
python 01_data_loading_and_merging.py
python 02_feature_cleaning_and_selection.py
python 03_model_training.py
```

Each accepts `--config <path>` to point at a non-default config. Paths in the
config are resolved relative to the folder that contains the `config/` directory.

What each script writes:

| Script | Main outputs |
|---|---|
| 01 | `data/TOPMed_merged.csv`, with columns in the order fixed by `config/feature_column_order.txt` |
| 02 | `data/TOPMed_cleaned.csv`, `data/final_feature_list.csv`, `data/TOPMed_gene_ids.csv`, `predict/training_medians.json` |
| 03 | `model/truncat.cbm` / `.pkl`, `results/` (CV predictions, importances, SHAP, figures) |

Script 01 needs the five annotation files in `data/annotations/`. They are
produced by the scripts in `Features/` and are not committed.

**Column order matters.** Results change by about 0.002 AUC with the column
order, so script 01 refuses to run if the merged columns differ from
`config/feature_column_order.txt`. If you change the feature set on purpose,
regenerate that file and rerun scripts 02 and 03.

## Check that the scripts reproduce the committed results

`check_parity.py` runs the scripts in a scratch directory (nothing in the repo
is overwritten) and compares their outputs with the committed files.

```bash
python check_parity.py                  # stages 01-02: fast, values must match exactly
python check_parity.py --through 03     # also trains the model (slow: cross-validation, final model, SHAP)
python check_parity.py --workdir /tmp/parity --no-run --through 03   # compare outputs of an earlier run
```

Stage 01 is skipped automatically if the annotation files are missing, and
stage 02 then starts from the committed `data/TOPMed_merged.csv`. For stage 03
the comparison uses tolerances (out-of-fold AUC within 0.003 by default).

For a full description of inputs, outputs, and operations at each stage,
see [`../notebooks/README.md`](../notebooks/README.md). The scripts mirror
the notebooks, and `check_parity.py` is how that is verified.

## Known differences from the notebooks

- The notebook `03_model_training.ipynb` also draws a feature-ablation curve
  (`feature_ablation_curve.png`). The script does not.
- The notebooks work from the notebook's own folder; the scripts resolve paths
  from the config location instead.
