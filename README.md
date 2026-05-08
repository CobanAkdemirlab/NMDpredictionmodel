NMDpredictionmodel/
├── 1_Data/
│   ├── TOPMed_extraction.R         
│   ├── gnomAD_extraction         
│   ├── ClinVar_extraction         
│   └── GREGoR_extraction         
│
├── Features/
│   ├── Canonical transcript annotation.R 
│   ├── gencodev26_features.R
│   ├── Variant annotation.R
│   ├── Transcript_Features.R
│   ├── simulation.R
│   ├── PTC_features.R
│   ├── motif_regions_extraction.py
│   ├── motif_fimo_matrix.py
│   ├── PTC_amino_acid.R
│   ├── Median_expression_mRNA Half-Life.R
│   ├── transcript_to_gene_mapping.py
│   ├── conservation_scores.py     # PhastCons & PhyloP
│   ├── ejc_occupancy_txnames_parallel_modified.py
│   ├── ptc_aug_analyzer.py
│   ├── conservation_scores.py     # PhastCons & PhyloP
│   ├── readthrough_scoring.py     # Readthrough potential
│   ├── ejc_analysis.py            # EJC occupancy
│   └── annotated_output/          # Annotated variant tables
│
├── 3_Model/

├── 4_Plotting/
│   ├── plot_predictions.py        # Prediction visualizations
│   ├── feature_distributions.py   # Feature plots
│   └── figures/                   # Output plots
│
├── 5_Interpretation/
│   ├── shap_analysis.py           # SHAP value generation
│   ├── shap_plots/                # SHAP visualizations
│   └── feature_importance/        # Importance rankings
└── README.md
