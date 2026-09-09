table truncat_predictions
"TrunCat NMD-escape predictions"
(
string  chrom;          "Chromosome"
uint    chromStart;     "Start position (0-based)"
uint    chromEnd;       "End position"
string  name;           "Variant ID (chrom_pos_ref_alt)"
uint    score;          "Escape probability x1000, for track shading (0-1000)"
char[1] strand;         "Strand (not available at variant level; '.')"
string  geneId;         "Ensembl gene ID"
string  hgncSymbol;     "HGNC gene symbol"
string  transcriptId;   "Transcript(s) used for prediction"
float   escapeProbability; "Raw predicted NMD-escape probability (0-1)"
string  predictedLabel; "Predicted class label (NMD / escape)"
float   thresholdUsed;  "Decision threshold used to assign predictedLabel"
string  cohort;         "Source cohort (gnomAD / ClinVar / GREGoR / TOPMed)"
string  model;          "Model used to generate this prediction"
)
