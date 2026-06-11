# NMDetective-B reimplementation

Applies the NMDetective-B decision tree (Lindeboom et al., Nat Genet
2019, Fig. 1c) to our stop-gain variant tables using locally-computed
features. This produces a head-to-head comparison against
TrunCat / TrunKitten on the same TOPMed held-out set, intended for the
manuscript supplement.

## Why reimplement instead of looking up figshare scores

The figshare resource (DOI 10.6084/m9.figshare.7803398) provides
per-position NMDetective scores annotated against UCSC knownGene
transcripts. Our pipeline annotates against GENCODE v26 (ENST IDs), so
direct lookup against the figshare GTFs falls back to any-overlap on
nearly every variant, defeating the point of a transcript-matched
comparison.

Reimplementing the published decision tree locally gives us:
- The same transcript context TrunCat saw — apples-to-apples evaluation
- No GENCODE↔UCSC mapping layer to maintain or defend in the methods
- Full transparency about the four rules being applied

NMDetective-A (Random Forest) is not reimplemented — it adds only ~3% R²
over NMDetective-B (71% vs 68% on the held-out indel set per Lindeboom
2019), and the paper itself uses NMDetective-B as the default in
downstream analyses, so this is the more relevant comparison anyway.

## The decision tree (Fig. 1c)

Evaluated in order; first match wins:

| # | Test                                       | Score | Rule label       |
|---|--------------------------------------------|-------|------------------|
| i | PTC in last exon                           | 0.00  | `last_exon`      |
| ii| Distance from PTC to coding start <150 nt  | 0.12  | `start_proximal` |
| iii| PTC exon length >407 nt                   | 0.41  | `long_exon`      |
| iv| PTC in last 50 nt of penultimate exon      | 0.20  | `50nt_rule`      |
| — | otherwise                                  | 0.65  | `trigger_NMD`    |

Higher score = more decay = more NMD-sensitive. Class assignment uses
the paper's published thresholds (Results section):

| Score          | Class          |
|----------------|----------------|
| > 0.52         | `sensitive`    |
| 0.25–0.52      | `intermediate` |
| < 0.25         | `escape`       |

## Required features

The script auto-detects column names from a candidate list, so it works
on `TopMed_merged_v4.csv` directly. Required fields:

| NMDetective feature                     | Default column         |
|-----------------------------------------|------------------------|
| on_last_exon                            | `last.exon`            |
| distance_to_coding_start (nt)           | `PTC.2.start`          |
| exon_length (nt)                        | `current.exon.length`  |
| in_last_50nt_of_penultimate_exon        | `penultimate.last50bp` |

Booleans tolerate R-style `TRUE`/`FALSE`/`T`/`F` as well as Python
`True`/`False` and 0/1.

## Usage

```bash
python nmdetective_b.py \
    --variants Model/TrunCat/data/TOPMed_merged.csv \
    --out      outputs/topmed_with_nmdetective_b.tsv \
    --qc       outputs/topmed_nmdetective_b_qc.tsv
```

For ~6,000 variants this runs in under a second. The output is the
input table plus three columns:

- `NMDetective_B_score` — continuous score from the leaf (0.00–0.65)
- `NMDetective_B_rule`  — which branch fired (one of the rule labels above)
- `NMDetective_B_class` — `sensitive` / `intermediate` / `escape` / `missing`

The QC tsv reports counts per rule and per class plus median/mean score.

## Tests

```bash
python test_nmdetective_b.py
```

Runs 9 cases covering every leaf, both strict-inequality boundary
conditions (150 nt and 407 nt thresholds), rule precedence (last_exon
dominates start_proximal), and missing-feature handling.

## Outputs

| Path                                       | Notes                  |
|--------------------------------------------|------------------------|
| `outputs/topmed_with_nmdetective_b.tsv`    | gitignored             |
| `outputs/topmed_nmdetective_b_qc.tsv`      | small; OK to commit    |

## Comparison against TrunCat / TrunKitten

Lives in a notebook (suggested:
`Model/TrunCat/notebooks/05_nmdetective_comparison.ipynb`), not here.

For the AUC comparison: TrunCat predicts probability of *escape*, while
NMDetective-B predicts NMD efficacy (higher = more decay). To score
NMDetective-B on your binary escape labels:

```python
from sklearn.metrics import roc_auc_score
# Lower NMDetective_B_score => more likely to escape
auc_nmd_b = roc_auc_score(y_escape, -df["NMDetective_B_score"])
```
