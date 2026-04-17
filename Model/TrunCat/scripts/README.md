# TrunCat — Standalone Scripts

Python script equivalents of the three notebooks. Same inputs, same outputs,
same config. Use when you want command-line execution instead of Jupyter.

## Requirements

- Python 3.12+
- Dependencies in `requirements.txt`

```bash
python -m venv venv
source venv/bin/activate       # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Verify setup

```bash
python setup.py
```

Checks config validity, input file presence, output directory creation, and
dependency installation.

## Usage

Run in order from `Model/TrunCat/scripts/`:

```bash
python 01_data_loading_and_merging.py
python 02_feature_cleaning_and_selection.py
python 03_model_training.py
```

Each accepts `--config <path>` to point at a non-default config.

For a full description of inputs, outputs, and operations at each stage,
see [`../notebooks/README.md`](../notebooks/README.md) — the scripts mirror
the notebooks exactly.