## External Validation Using Saturation Genome Editing (SGE) Data

To evaluate the applicability of TrunKitten to experimentally characterized variants, we applied the trained model to an independent saturation genome editing (SGE) dataset for *LMNA*.

The SGE variants were processed to generate the 10 features required by the trained TrunKitten model. The annotated feature matrix was then used to generate predictions without retraining the model.

### Running TrunKitten Predictions

Run the following command from the repository root:

```bash
python Model/TrunKitten/pipeline/predict.py \
  --annotated output/LMNA_TrunKitten_10features_with_id.tsv \
  --model Model/TrunKitten/model/trunkitten.pkl \
  --metadata Model/TrunKitten/model/trunkitten_features.json \
  --out output/LMNA_TrunKitten_predictions.tsv
```

### Input and Model Files

| File | Description |
|---|---|
| `LMNA_TrunKitten_10features_with_id.tsv` | Annotated LMNA SGE variants containing the 10 model-required features and variant identifiers. |
| `trunkitten.pkl` | Pretrained TrunKitten prediction model. |
| `trunkitten_features.json` | Model feature metadata used to ensure compatibility between the input data and the trained model. |

### Output

The prediction results are saved to:

`output/LMNA_TrunKitten_predictions.tsv`

This file contains the TrunKitten predictions for the annotated LMNA SGE variants.

### External Validation

The predictions can be compared with experimentally measured SGE outcomes to assess concordance between TrunKitten predictions and observed variant effects.

The SGE dataset is used for independent evaluation and is not included in model training.
