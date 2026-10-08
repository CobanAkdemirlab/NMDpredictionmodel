# TrunCat UCSC track hub

A UCSC Genome Browser track (hg38) of TrunCat NMD-escape predictions for the
scored variants in gnomAD, ClinVar and GREGoR. It covers scored variants only,
not the whole genome.

**Load the track in the UCSC Genome Browser:**
[TrunCat NMD predictions (hg38)](https://genome.ucsc.edu/cgi-bin/hgTracks?hubUrl=https://raw.githubusercontent.com/CobanAkdemirLab/NMDpredictionmodel/main/Model/TrunCat/track/hub.txt&db=hg38)

Or add it manually: in the Genome Browser go to *My Data -> Track Hubs*, open
the *Connected Hubs* tab, and paste

```text
https://raw.githubusercontent.com/CobanAkdemirLab/NMDpredictionmodel/main/Model/TrunCat/track/hub.txt
```

## Files

| File | Purpose |
|------|---------|
| `hub.txt`, `genomes.txt`, `trackDb.txt` | Hub definition read by UCSC |
| `truncat_predictions.bb` | The track (bigBed) |
| `truncat_predictions.bed`, `truncat_predictions.as` | Source BED and its field definitions |
| `hg38.chrom.sizes` | Chromosome sizes used to build the bigBed |

## Fields

Each feature is one scored variant. Fields include the variant ID
(`chrom_pos_ref_alt`), Ensembl gene ID and gene symbol, transcript(s), the
escape probability (0-1), a predicted class (`NMD` or `escape`), the decision
threshold used, and the source cohort. Feature shading is the escape
probability x 1000. The complete field list is in `truncat_predictions.as`.

A variant that appears in more than one cohort has one feature per cohort.

## Rebuilding

The track is built from the cohort prediction files in `../predict/` by
`../notebooks/build_ucsc_track.ipynb`. Rebuild and recommit the `.bed` and
`.bb` files whenever the predictions change.
