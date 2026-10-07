# Pipeline dependencies

This directory contains dependencies required for the upstream
NMDpredictionmodel workflow:

1. dataset extraction;
2. variant and transcript annotation;
3. biological feature generation.

These dependencies are separate from the machine-learning dependencies
used for TrunCat model development.

## Python

Install upstream Python dependencies with:

```bash
pip install -r dependencies/requirements_pipeline.txt
