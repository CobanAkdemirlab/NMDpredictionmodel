# Features

This directory contains the scripts used to generate genomic, transcriptomic, sequence-context, and regulatory features for TrunCat. The feature-generation framework is shared across TOPMed, gnomAD, ClinVar, and GREGoR after dataset extraction and transcript/PTC annotation.

**Scope:** This README describes the purpose and major dependencies of each feature family. Exact input/output filenames and command-line arguments should be taken from the corresponding script's configuration or help text; the scripts do not all use the same filename conventions.

## Workflow

```text
Dataset_extraction/ + Annotation/
                 |
                 v
      Annotated PTC variants
                 |
                 v
     Feature-family scripts (01-13)
                 |
                 v
     Feature tables / motif matrix
                 |
                 v
       99_build_feature_matrix.R
                 |
                 v
       Model/TrunCat/ prediction
```

Feature families are generally independent; numbered filenames organize them, but do not imply that every feature script must be run strictly in numeric order. Some families have required internal substeps (for example, `07a` before `07b`).

## Feature-generation scripts

| Script | Feature family | Main inputs or references | Purpose |
|---|---|---|---|
| `01a_GENCODEv26_sequence_features.R` | Transcript sequence | Annotated variants, GENCODE v26 transcript/CDS/UTR sequences | Sequence lengths and nucleotide composition |
| `01b_GENCODEv26_intron_features.py` | Transcript introns | GENCODE v26 annotation, transcript identifiers | CDS, 5′UTR, and 3′UTR intron annotations |
| `02_PTBP1_binding_features.R` | PTBP1 binding | Transcript regions and PTBP1 binding intervals | Binding-related features in transcript regions |
| `03_PTC_amino_acid_context.R` | PTC context | PTC and coding-sequence annotations | Amino-acid context around the PTC |
| `04_gene_level_features.R` | Gene properties | Gene identifiers, GTEx expression, gnomAD constraint, mRNA stability data | Gene-level expression, constraint and half-life covariates |
| `05_PTC_geometry_derived_features.R` | PTC geometry | PTC position, CDS/exon structure | Relative PTC position, exon-boundary distances and related NMD geometry |
| `06a_prepare_variant_scores.R` | Variant scores: preparation | Genomic variant identifiers | Prepare external variant-score annotation input |
| `06b_run_variant_scores.sh` | Variant scores: annotation | Prepared variants and configured ANNOVAR resources | Run score annotation |
| `06c_merge_variant_scores.R` | Variant scores: merge | ANNOVAR score outputs | Merge variant-level scores with variant identifiers |
| `07a_motif_region_extraction.py` | Motif regions | PTC/transcript annotations and sequence reference | Extract regions around PTCs, exon junctions and UTRs |
| `07b_FIMO_motif_features.py` | RBP motifs | Extracted sequences, ATtRACT motifs, MEME Suite/FIMO | Scan sequences and build motif presence features |
| `08a_make_optimal_codons_from_trna.py` | Codon reference | hg38 UCSC/GtRNAdb-derived tRNA table | Construct optimal-codon reference list |
| `08b_codon_optimality_features.py` | Codon optimality | Coding sequences and optimal-codon list | CDS-wide and local codon-optimality features |
| `09_conservation_score_features.py` | Conservation | hg38 phastCons100way and phyloP100way BigWigs | Conservation scores in relevant sequence regions |
| `10_EJC_occupancy_features.py` | EJC occupancy | EJC occupancy reference intervals (GSE41154) | Occupancy-related features near exon junctions |
| `11_PTC_AUG_features.py` | Downstream AUG / Kozak | Transcript sequences and PTC coordinates | Downstream AUG frames, distances and Kozak context |
| `12_readthrough_features.py` | Readthrough | PTC nucleotide sequence context | Stop-codon readthrough-context features |
| `13_isoforn_count_features.py` | Isoform counts | Transcript/gene annotations | Transcript-isoform count features |
| `99_build_feature_matrix.R` | Feature assembly | Annotated variants and generated feature tables | Assemble model-ready feature data |

The script name `13_isoforn_count_features.py` is reproduced as listed in the repository documentation; confirm the filename in GitHub if it has been renamed.

## Dependencies between scripts

These subworkflows have a specific internal order:

```text
06a_prepare_variant_scores.R
          |
          v
06b_run_variant_scores.sh
          |
          v
06c_merge_variant_scores.R

07a_motif_region_extraction.py
          |
          v
07b_FIMO_motif_features.py

08a_make_optimal_codons_from_trna.py
          |
          v
08b_codon_optimality_features.py
```

Other feature scripts require their own annotation/reference inputs and may be run independently after those inputs are prepared.

## Reference resources

| Resource | Used for | Notes |
|---|---|---|
| GENCODE v26 GTF and transcript sequences (GRCh38) | Sequence, intron, PTC and transcript features | Transcript identifiers and versions should be consistent across stages |
| GTEx v8 expression | Gene-level features | Configure reference path in the gene-level script |
| gnomAD v4.1 constraint metrics | pLI / loss-of-function constraint | Canonical gene rows were used to avoid duplicated gene joins |
| mRNA half-life data / PC1 | Gene-level mRNA stability | External derived reference; document its location |
| ATtRACT motif definitions | RNA-binding protein motif features | Used with FIMO; the analysis used a p-value cutoff of 1e-4 and sense-only scanning |
| `reference/codon_optimality/hg38_UCSC_tRNA_table.tsv` | Optimal-codon derivation | Input to `08a` |
| `reference/codon_optimality/optimal_codons.txt` | Codon-optimality features | Output/reference consumed by `08b` |
| hg38 phastCons100way and phyloP100way BigWigs | Conservation | Large reference files may need separate downloads |
| GEO GSE41154 | EJC occupancy | Confirm genome-build compatibility before using interval coordinates |
| ANNOVAR annotation databases | Variant scores | Requires a locally configured ANNOVAR installation and resources |

Some large or licensed references may not be distributed with the repository. The presence of a script does not guarantee its reference file is included in GitHub.

## Feature matrix and prediction

`99_build_feature_matrix.R` is intended to integrate feature families into a common variant-level table. The TrunCat classifier uses **619 selected model features** from the TOPMed model-development analysis (**5,749 variants**).

For external prediction, align input columns to the **trained model's feature names and types** using the prediction utilities in [`../Model/TrunCat/predict/`](../Model/TrunCat/predict/). Missing-value handling should follow the trained preprocessing workflow rather than ad hoc substitutions. ClinVar `CLNSIG` and `clinsig_group` are retained as metadata and must not be used as model predictors.

## Running and reproducibility notes

1. Complete [`Dataset_extraction/`](../Dataset_extraction/README.md) and [`Annotation/`](../Annotation/README.md) first.
2. Configure genome/transcript resources and script-specific input/output directories.
3. Run the needed feature families, respecting the dependent substeps above.
4. Check variant keys, transcript identifiers, row counts, missing values, and generated columns before feature assembly.
5. Assemble and align features for the trained TrunCat model, then predict using [`Model/TrunCat/`](../Model/TrunCat/).

**Important:** The feature scripts and matrix assembly have had compatibility issues in earlier revisions (including join keys and expected output filenames). The repository-wide workflow has not been independently verified as a single clean end-to-end run. Users should validate outputs at each stage; do not interpret this README as a claim that all scripts currently run without modification.
