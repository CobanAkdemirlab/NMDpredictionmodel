# TrunKitten Annotation Pipeline — Design & Specification

**TrunKitten** is the reduced top-10 feature NMD-prediction model derived
from **TrunCat** (TRUNcation-aware Classifier using Annotated Transcripts).
This pipeline annotates pre-called PTC / stop-gain variants with exactly
the 10 features TrunKitten requires. All conventions are aligned with the
TrunCat training-pipeline scripts (`conservation_score_extraction_v2.py`,
`ejc_occupancy_txnames_parallel_modified.py`, and feature-generation code
in the `CobanAkdemirLab/NMDpredictionmodel` repo), so TrunKitten predictions
on external cohorts are directly comparable to the TrunCat training cohort.

> The implementation package is named `minicat` for historical reasons;
> all module paths and imports remain unchanged.

**Do not use this pipeline to call PTCs.** It assumes the input rows are
already confirmed stop-gained SNVs with a valid per-variant `txnames`
transcript assignment.

---

## 1. Feature definitions

All transcript-oriented coordinates use **1-based, inclusive** positions within
the spliced mature transcript (5′ UTR + CDS + 3′ UTR), ordered 5′→3′ in the
transcript's reading direction. Genomic coordinates are 1-based inclusive on
input (GTF convention) and converted to 0-based half-open only at the
BED/BigWig boundary.

### Universal convention: the "PTC position"

The variant `position` column is the **genomic coordinate of the mutated base**
(the SNV ref→alt site). Everywhere downstream, "PTC" = this genomic position.
This matches `conservation_score_extraction_v2.py` exactly.

**Ambiguity: "after the PTC" = after the full stop codon or after the first mutated base?**
Plausible interpretations:
- (a) After the first mutated base (genomic `position`).
- (b) After the full 3-nt stop codon (position + 2 on + strand, position − 2 on − strand).

**Recommendation: (a), the first mutated base.** This matches the existing
training-data conservation script (`new3utr_blocks` uses `ptc + 1` on +strand
and `ptc − 1` on −strand — one base past the mutated site, which is the SNV
position). Using (b) would shift all `new3utr` regions by 2 bp relative to the
training distribution and create silent skew in `phastcons_new3utr_first200_median`.

### 1.1 `last.EJC` — categorical {`upstream`, `penultimate.last50bp`, `last.exon`}

- **Meaning**: Canonical NMD positional rule. PTCs in the last exon or within
  50 nt of the last exon-exon junction typically escape the canonical EJC-EJC
  mechanism.
- **Rule**: Use **all transcript exons** (not coding-only), ordered 5′→3′ along
  the transcript:
  - `last.exon` ⇢ PTC is in the final exon (highest transcript-order rank).
  - `penultimate.last50bp` ⇢ PTC is in the second-to-last exon **AND** the
    transcript-oriented distance from the PTC to the 3′ end of that exon is
    ≤ 50 nt.
  - `upstream` ⇢ everything else.
- **Exon ordering**: On + strand, 5′→3′ = ascending genomic start. On −
  strand, reverse — the exon with the smallest genomic start is the *last*
  exon. We compute rank from strand-oriented order.
- **Last 50 bp definition (penultimate exon, transcript-oriented)**:
  - + strand: `penultimate.end - ptc ≤ 50` (genomic `end` is the 3′ boundary).
  - − strand: `ptc - penultimate.start ≤ 50` (genomic `start` is the 3′ boundary).
- **Boundary PTCs**: If the PTC genomic position sits exactly at an exon
  boundary, it is assigned to the exon whose interval includes it (exon
  intervals are closed `[start, end]`). A variant at a splice site itself is
  flagged in the QC `boundary_ambiguous` column but still assigned by the
  inclusive rule above.
- **Output**: string in {`upstream`, `penultimate.last50bp`, `last.exon`}.
- **Edge cases**: Single-exon transcripts → always `last.exon`.

### 1.2 `mut.exon` — integer

- **Meaning**: The 1-based transcript-order rank of the PTC-containing exon,
  across **all** transcript exons (not coding-only). This is the standard
  "exon number" that matches the training feature (redundant with
  `AmountExonsBefore`, which was dropped because `mut.exon` encodes the same
  information).
- **Ordering**: strand-oriented 5′→3′, same as §1.1.
- **Boundary handling**: Same inclusive rule as `last.EJC`.
- **Output**: int ≥ 1.

### 1.3 `relativePTClocation` — float in [0, 1]

- **Numerator**: transcript-oriented nucleotide position of the PTC within the
  full spliced mature transcript (1-based). Computed by walking exons
  5′→3′ and summing lengths of upstream exons plus the within-exon offset.
- **Denominator**: total spliced transcript length = Σ(exon lengths) (UTR +
  CDS + UTR; matches the training denominator).
- **Formula** (+ strand):
  ```
  transcript_pos = Σ(len of exons k where k < mut.exon) + (ptc - exon_start + 1)
  relativePTClocation = transcript_pos / Σ(all exon lengths)
  ```
- **Formula** (− strand):
  ```
  transcript_pos = Σ(len of exons k where k < mut.exon) + (exon_end - ptc + 1)
  ```
- **PTC convention**: uses genomic `position` (first mutated base); same
  convention as `conservation_score_extraction_v2.py`.
- **UTR included**: yes — denominator is total exonic length.
- **Output**: float in (0, 1].

### 1.4 `AmountExonsAfter` — integer ≥ 0

- **Meaning**: number of **coding exons** strictly downstream (in transcript
  order) of the exon that contains the PTC.
- **Coding exon identification**: A transcript exon is "coding" if it
  overlaps any CDS record for that transcript in the GTF, i.e.
  `exon.start ≤ cds.end` AND `exon.end ≥ cds.start` for any CDS entry. Use
  Gencode `feature=="CDS"` records.
- **PTC-containing exon**: counted as the reference exon; not included in the
  "after" count even if only partially coding (e.g. the PTC sits in a
  partially-coding exon — it is the *current* exon, not an after-exon).
- **Transcript order**: strand-oriented 5′→3′, same as §1.1.
- **Last-exon PTCs**: `AmountExonsAfter = 0`.
- **Output**: int ≥ 0.

### 1.5 `cdsseqs_AU_content` — float in [0, 1]

**"cdsseqs" — disambiguation.**

Plausible interpretations:
- (a) The **full** annotated CDS sequence of the transcript (start codon → native stop).
- (b) The **truncated** CDS from start codon up to and including the PTC codon.
- (c) The truncated CDS up to but excluding the PTC codon.

**Recommendation: (a), the full annotated CDS**, start codon through (and
including) the native stop codon. Rationale:
1. The training feature name is plural "cdsseq**s**" — plural = the full CDS
   strings as annotated, not a variant-specific truncation.
2. CDS-level composition features in the training set (`cdsseqs_AU_content`,
   `cdsseqs_UC_content`) correlate with transcript-intrinsic half-life /
   codon-optimality signals, which are properties of the annotated mRNA —
   not of the mutated, truncated product.
3. This makes the feature depend only on `txnames`, not on the variant
   position, so it's stable across any variant hitting the same transcript.

**DNA vs RNA alphabet.** Reference FASTA stores DNA (A/C/G/T). For
transcript-oriented sequence on − strand we reverse-complement the genomic
sequence. AU/UC content is a *biological* (RNA) concept but equivalent to
AT / TC on the DNA alphabet. We compute on the DNA alphabet and name the
feature as-trained (AU, UC) — no character substitution needed:

- `AU_content = (count(A) + count(T)) / total_non_N_bases`
- `UC_content = (count(T) + count(C)) / total_non_N_bases`

"Content" = fraction of single bases in the set, **not** dinucleotides.
N bases are excluded from numerator and denominator.

### 1.6 `cdsseq_AUcontentlast200` — float in [0, 1]

Note: singular `cdsseq` here vs plural `cdsseqs` above. The existing feature
files in the repo use this naming distinction; we preserve it.

**"last200" disambiguation.**

Plausible interpretations:
- (a) Last 200 nt of the annotated CDS (before the native stop codon).
- (b) Last 200 nt before the PTC (i.e. the last 200 nt of the truncated CDS).
- (c) Last 200 nt of the full CDS including the native stop codon.

**Recommendation: (a), the last 200 nt of the annotated CDS immediately
before the native stop codon.** Rationale:
1. The zero-fill logic in `02_feature_cleaning_and_selection.ipynb`
   zero-fills this feature only when the CDS is shorter than the window size
   — confirming it's a fixed annotation-level window, not a variant-truncated
   one. A variant-truncated "last 200 before PTC" would be missing for any
   variant too close to the CDS start, which is *not* how the training
   feature is distributed.
2. This feature appears in rows of the conservation pipeline keyed by
   transcript only; its value would be identical for two PTCs in the same
   transcript — consistent with (a).
3. Interpretation (a) makes the feature independent of variant position,
   matching the behavior of `cdsseqs_AU_content`.

**Definition:** Take the transcript-oriented CDS sequence (start codon →
base just before native stop codon, **excluding** the stop codon itself).
Take the last 200 bases of this sequence (or the full CDS if `len < 200`, and
emit `cds_length_short_flag=True`). Compute `AU_content` on those bases.

If CDS length < 200, the training pipeline zero-fills this (see
`utr_zero_fill` in Notebook 02). We return the actual fraction and flag it;
the inference user applies the same zero-fill rule before scoring.

### 1.7 `cdsseqs_UC_content` — float in [0, 1]

Same as §1.5 but `count(T) + count(C)` / `total_non_N`. "UC" in the RNA
alphabet = T+C on DNA. Computed on the full annotated CDS (§1.5 convention).

### 1.8 `half_life_PC1` — float

- **Source**: Excel file indexed by ENSG gene IDs.
- **Merge logic**:
  1. From the GTF, for each `txnames` (transcript ID), read the transcript
     feature's `gene_id` attribute.
  2. Strip the Gencode version suffix (`ENSG00000123456.5` → `ENSG00000123456`).
  3. Strip version suffixes from the Excel table's ENSG keys too.
  4. Left-join.
- **Missing handling**: retain `NaN`. The training pipeline median-imputes
  `half_life_PC1` (median = 1.1303 in the v4 cleaned file; see Notebook 02
  log). For external prediction we output `NaN` and a `half_life_PC1_missing`
  flag; the inference-time consumer applies the same median from the
  training distribution. (Imputing per-cohort medians would be wrong — it
  would drift the feature distribution relative to training.)

### 1.9 `phastcons_new3utr_first200_median` — float

- **Region** (transcript-oriented): from the base immediately **after** the
  PTC position through the 3′ end of the transcript, truncated to the first
  200 transcript bases.
- **"After the PTC" convention**: first mutated base + 1 (see universal
  convention above). This matches `new3utr_blocks` in
  `conservation_score_extraction_v2.py`.
- **Mapping to genomic coordinates**: walk transcript exon blocks 5′→3′,
  skip blocks entirely upstream of the PTC, clip the PTC-containing block
  to `(ptc+1 .. exon_end)` on + strand or `(exon_start .. ptc-1)` on −
  strand, then take the first 200 transcript bases across subsequent exon
  blocks (uses `take_first_N` logic from the existing script — strand-aware
  block ordering + last-block clipping).
- **Median computation**: concatenate BigWig values across all constituent
  genomic intervals, drop NaN, take `np.median`. If every base is missing,
  return `NaN`.
- **Last-exon PTCs**: the new3utr region exists (it spans the remainder of
  the last exon). Generally shorter than 200 nt; `downstream_new3utr_len`
  records the actual length.
- **Missing / absent region**: `NaN`, plus `new3utr_empty` flag. Training
  pipeline zero-fills; inference consumer applies the same zero-fill.
- **BigWig vs BED input**:
  - BigWig → `pyBigWig.values(chrom, start0, end0, numpy=True)` per segment.
  - BED/bedGraph → build an interval tree (e.g. via `pyranges`) keyed by
    chromosome; for each region segment, compute the weighted per-base
    vector by intersecting with BED entries and carrying the score field.
    We provide a `_extract_from_bed` fallback in `conservation.py`.

### 1.10 `phylop_ptc_to_ejc_median` — float

- **Region** (genomic, within the PTC-containing exon, strand-aware):
  - + strand: `[ptc+1 .. exon_end]` of the current exon.
  - − strand: `[exon_start .. ptc-1]` of the current exon.
  - This is the **within-exon** downstream stretch up to the nearest
    exon-exon junction. Matches `conservation_score_extraction_v2.py`
    exactly.
- **"Downstream EJC" definition**: the nearest 3′ exon-exon junction in
  transcript orientation, which is the current exon's 3′ boundary. The
  region only covers up to that junction — not across it.
- **Special cases**:
  - `last.EJC == 'last.exon'` → no downstream EJC exists → region is
    empty → `NaN` (training pipeline zero-fills).
  - `last.EJC == 'penultimate.last50bp'` → region still exists (penultimate
    exon has a downstream EJC at its 3′ end). The "last 50 bp" rule is
    about NMD biology, not about whether the region exists.
- **Missing / no-valid-interval**: `NaN`, plus `ptc_to_ejc_empty` flag.
- **Implementation**: identical median-over-valid-values logic as §1.9.

---

## 2. Pipeline architecture

### Language & dependencies

Python ≥ 3.10. Install into a fresh venv:

```bash
pip install pandas numpy pyfaidx pyBigWig pyranges gffutils openpyxl catboost joblib
```

Optional: `pysam` (not strictly needed — `pyfaidx` handles FASTA; `pyranges`
handles BED).

### Directory layout

```
minicat_pipeline/
├── config/
│   └── config.yaml              # paths, conventions, flags
├── inputs/                      # user-supplied
│   ├── variants.tsv
│   ├── annotation.gtf.gz
│   ├── genome.fa(.fai)
│   ├── phastcons.bw
│   ├── phylop.bw
│   └── half_life_pc1.xlsx
├── minicat/                     # package
│   ├── __init__.py
│   ├── config.py                # yaml loader
│   ├── gtf_index.py             # build exon/CDS index by transcript
│   ├── transcript.py            # position mapping, exon ranking
│   ├── sequence.py              # CDS reconstruction, AU/UC content
│   ├── conservation.py          # BigWig / BED median extraction
│   ├── halflife.py              # Excel merge by ENSG
│   ├── features.py              # per-variant feature assembly
│   ├── qc.py                    # QC columns + missing flags
│   └── cli.py                   # CLI entry point
├── outputs/
│   ├── annotated.tsv            # final feature table
│   ├── qc_report.tsv            # per-variant QC
│   └── run.log
└── tests/
    └── test_toy_transcripts.py  # §5 validation
```

### Module responsibilities

| Module           | Responsibility                                                                 |
| ---------------- | ------------------------------------------------------------------------------ |
| `config.py`      | Load YAML, validate paths exist, return immutable config dataclass.            |
| `gtf_index.py`   | Parse GTF once → `{transcript_id: TranscriptRecord}` with exons, CDS, gene_id. Normalises version suffixes.  |
| `transcript.py`  | `genomic_to_transcript_pos`, `locate_exon`, `rank_exons_5to3`, `coding_exons`. |
| `sequence.py`    | Build spliced transcript string (strand-aware RC), extract CDS, compute AU/UC. |
| `conservation.py`| Extract `new3utr_first200` and `ptc_to_ejc` intervals, query BigWig/BED, return median + QC length. |
| `halflife.py`    | Load Excel; map `txnames → ENSG → half_life_PC1`.                              |
| `features.py`    | Orchestrator: for each variant row, produce the 10-feature dict + QC fields.   |
| `qc.py`          | Assembles missing-flag columns, boundary-ambiguous markers, length QC.         |
| `cli.py`         | Argparse entry point; parallelises per-variant annotation; writes outputs.     |

### Logging / QC strategy

- Python `logging` → both stderr and `outputs/run.log` at INFO.
- Per-variant failures are collected into a warnings list, not raised.
- `qc_report.tsv` has one row per input variant with columns:
  `variant_id, txnames, transcript_id_used, gene_id, strand, exon_count,
  coding_exon_count, ptc_transcript_pos, transcript_length, cds_length,
  downstream_new3utr_len, boundary_ambiguous, tx_not_in_gtf, version_mismatch,
  cds_length_short_flag, new3utr_empty, ptc_to_ejc_empty, half_life_missing,
  any_conservation_missing`.
- Final summary at exit: count of successfully annotated variants, missing
  transcripts, missing half-life keys, empty regions.

---

## 3. Region / sequence logic (consolidated)

| Region                | Transcript-oriented definition                              | Genomic realisation                                                                                       |
| --------------------- | ----------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| new3utr_first200      | First 200 tx-bases strictly after the PTC position          | Clip PTC-containing exon at (ptc±1, boundary); concat downstream exon blocks in transcript order; take_first_N=200. |
| ptc_to_ejc            | PTC+1 → nearest downstream exon-exon junction (within exon) | Single genomic interval within the PTC-containing exon: `[ptc+1, exon_end]` (+) or `[exon_start, ptc-1]` (−). |
| Full CDS (for content)| start codon → base before stop codon                        | Concatenate CDS GTF intervals in transcript order; RC if − strand.                                        |
| CDS last 200 (content)| Last 200 tx-bases of full CDS                               | Slice `[-200:]` of the reconstructed CDS string.                                                           |

BigWig mechanics: all regions are converted to BED half-open `[start0, end0)`
segments only at the query boundary. Medians are computed over the
concatenation of `bw.values()` arrays across all segments for a given region,
dropping NaN. No re-weighting by segment length (per-base values are
intrinsically per-base).

BED fallback: if the conservation file is a bedGraph-like BED with a `value`
column, we build per-chromosome `pyranges` and do an interval intersection,
expanding to per-base values using the entry's score. Median over valid bases.
(For phyloP / phastCons, BigWig is strongly preferred — bedGraph is order-of-
magnitude slower for per-base median.)

---

## 4. Starter code — see `minicat/` package

Each module is in `minicat/*.py`. See code for docstrings and examples.

---

## 5. Validation plan

Three toy transcripts are constructed programmatically in
`tests/test_toy_transcripts.py`:

### Toy 1 — `+` strand, 3 exons

```
      exon1: 100-199   exon2: 300-399   exon3: 500-699 (last)
      |--------|        |--------|       |---------------|
 CDS:         150-199,  300-399,         500-599    stop at 600-602
 tx strand = +
```

Variant A: `position = 350` (in exon 2).
- `last.EJC` = `upstream` (exon 2 is not last, and is not penultimate in a
  3-exon transcript — wait, it *is* the penultimate; distance from PTC to
  `end` = 399 − 350 = 49 → **`penultimate.last50bp`**). Good boundary test.
- `mut.exon` = 2.
- `AmountExonsAfter` = 1 (exon 3 is coding: 500-599 overlaps CDS).
- `relativePTClocation` = (100 + (350−300+1)) / (100+100+200) = 151/400 = 0.3775.

Variant B: `position = 340`.
- Distance to exon2 end = 59 → `upstream`.

Variant C: `position = 550` (in exon 3).
- `last.EJC` = `last.exon`; `mut.exon` = 3; `AmountExonsAfter` = 0.

### Toy 2 — `−` strand, 3 exons, same coordinates

Exon order flips: genomic-small-start exon is the *last* exon. Validates that
rank assignment produces `last.exon` for a PTC at `position=150` (in the
genomically-leftmost exon, which is the *last* exon on − strand).

### Toy 3 — Single-exon transcript

Always `last.exon`, `AmountExonsAfter = 0`, `mut.exon = 1`.

### Conservation region validation

Mock BigWig: build a small in-memory array over a known region, write with
`pyBigWig.open(..., "w")`, validate that our median matches a manually
computed median.

### Strand correctness

For each toy transcript, compare our reconstructed CDS string against a
hand-computed one (write out exon sequences manually).

---

## 6. Common failure modes / QC checks

| Failure mode                                           | Detection                                       | Behaviour                                                                                    |
| ------------------------------------------------------ | ----------------------------------------------- | -------------------------------------------------------------------------------------------- |
| `txnames` not in GTF (even after version strip)        | `tx_not_in_gtf=True`                            | Skip variant; emit NaN row with flag; log warning.                                           |
| `txnames` version differs from GTF                     | Normalise both sides by stripping `.N`          | Report `version_mismatch=True` (informational only).                                         |
| Input `position` not inside any exon of `txnames`      | Exon-containing lookup returns empty            | Skip variant; emit NaN row; log. This is a data quality issue — PTC should be exonic.        |
| PTC at exon boundary                                   | `position == exon_start` or `position == exon_end` | `boundary_ambiguous=True`; still assigned by inclusive rule.                              |
| `cds_length < 200` for `cdsseq_AUcontentlast200`       | Measured CDS length < 200                       | Return actual fraction; set `cds_length_short_flag=True`. Consumer applies training's zero-fill. |
| new3utr region empty (shouldn't happen if not last exon — but sanity check) | `take_first_N` returns empty list     | Set `new3utr_empty=True`; features → NaN.                                                    |
| `ptc_to_ejc` empty (last.exon)                         | `last.EJC == 'last.exon'`                       | Expected; `ptc_to_ejc_empty=True`; feature → NaN.                                            |
| All BigWig values NaN in region                        | `valid_bp == 0`                                 | Feature → NaN; `any_conservation_missing=True`.                                              |
| Excel missing ENSG                                     | Merge miss                                      | `half_life_missing=True`; feature → NaN.                                                     |
| Chromosome naming mismatch (`chr1` vs `1`) between GTF and BigWig | Harmonisation step at BigWig open | Auto-detect by probing `bw.chroms()` and prefixing / stripping `chr`.                        |

### Final output schema

One row per input variant, column order:

1. `variant_id` *(str)*
2. `txnames` *(str)* — as provided
3. `transcript_id_used` *(str)* — the GTF-matched id (may differ by version suffix)
4. `gene` *(str)* — as provided
5. `gene_id` *(str)* — ENSG from GTF (version-stripped)
6. `strand` *(str)* — `+` or `−`
7. **Features (10):**
   - `last.EJC` *(str / category)*
   - `relativePTClocation` *(float64)*
   - `half_life_PC1` *(float64 or NaN)*
   - `cdsseqs_AU_content` *(float64)*
   - `mut.exon` *(int)*
   - `phastcons_new3utr_first200_median` *(float64 or NaN)*
   - `phylop_ptc_to_ejc_median` *(float64 or NaN)*
   - `AmountExonsAfter` *(int)*
   - `cdsseq_AUcontentlast200` *(float64)*
   - `cdsseqs_UC_content` *(float64)*
8. **QC** *(all written to `qc_report.tsv`)*: `exon_count`, `coding_exon_count`,
   `ptc_transcript_pos`, `transcript_length`, `cds_length`, `downstream_new3utr_len`,
   `boundary_ambiguous`, `tx_not_in_gtf`, `version_mismatch`, `cds_length_short_flag`,
   `new3utr_empty`, `ptc_to_ejc_empty`, `half_life_missing`, `any_conservation_missing`.

Example row (tab-separated, abbreviated):
```
variant_id           txnames             transcript_id_used   gene   gene_id           strand  last.EJC               relativePTClocation   half_life_PC1   cdsseqs_AU_content   mut.exon   phastcons_new3utr_first200_median   phylop_ptc_to_ejc_median   AmountExonsAfter   cdsseq_AUcontentlast200   cdsseqs_UC_content
chr8_41977233_C_A    ENST00000265713     ENST00000265713.8    KAT6A  ENSG00000083168   -       upstream               0.3421                1.0523          0.511                5          0.982                                1.845                      12                 0.495                     0.498
```

---
