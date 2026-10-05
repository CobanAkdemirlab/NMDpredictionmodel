# Figure 2 Reproduction Guide

The analyses and plotting code used to generate **Figure 2** are available in:

`Model/TrunCat/notebooks/03_model_training.ipynb` or `Model/TrunCat/notebooks/recreate_shap_summary_from_pkl.ipynb`

The table below indicates the notebook section corresponding to each panel.

| Figure panel | Analysis | Notebook section |
|---|---|---|
| **Figure 2A** | SHAP value distributions for the 12 most influential features  | **Model/TrunCat/notebooks/recreate_shap_summary_from_pkl.ipynb** |
| **Figure 2B** | Predicted probability distribution | **Model/TrunCat/notebooks/03_model_training.ipynb (4. Probability Distribution)** |
| **Figure 2C** | Confusion matrix and classification performance | **Model/TrunCat/notebooks/03_model_training.ipynb (8. Confusion Matrix and Classification Report)** |
| **Figure 2D** | Top 20 predictors ranked by mean absolute SHAP value | **Model/TrunCat/notebooks/03_model_training.ipynb (7a. Top 20 Feature Importances by Mean \|SHAP\| (Primary))** |
| **Figure 2E** | Feature ablation analysis | **Model/TrunCat/notebooks/03_model_training.ipynb (9. Ablation Curve Report)** |
| **Figure 2F** | Comparison of TrunCat, TurnKitten and NMDetective-B ROC and precision-recall (PR) performance | **Model/Benchmarking/NMDetective-B/benchmark_nmdetectiveB.ipynb** |

## Notes

For **Figure 2C**, Section 8 generates the confusion matrix and classification metrics. These results are used to produce the summarized classification-performance panel shown in the manuscript.

For **Figure 2D**, feature importance is based on the **mean absolute SHAP value (mean |SHAP|)** across model predictions.

The final Figure 2 panels generated for the manuscript are also provided under `Model/TrunCat/` for reference.
