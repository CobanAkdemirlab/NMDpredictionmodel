# Reference inputs for the TrunKitten annotation pipeline

Paths are set in `../pipeline/config/config.yaml`. Files marked "in repo" are committed;
everything else is a public download that you place here (or point the config at
wherever you keep it). The large files are git-ignored.

| File (default name) | Source | In repo? |
| --- | --- | --- |
| `annotation.gtf.gz` | GENCODE **v26** (hg38) comprehensive/primary-assembly annotation, <https://www.gencodegenes.org/human/release_26.html>. Use v26: it matches TrunCat training and the GTEx v8 gene IDs. | no |
| `genome.fa` (+ `.fai`) | hg38 primary assembly FASTA, indexed with `samtools faidx genome.fa` | no |
| `phastcons.bw` | UCSC hg38 phastCons 100-way vertebrate bigWig (`hg38.phastCons100way.bw`), <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/> | no |
| `half_life_pc1.xlsx` | mRNA half-life PC1 from Agarwal & Kelley (2022), *Genome Biol* 23:245 (supplementary file `13059_2022_2811_MOESM3_ESM.xlsx`). The committed copy already has the columns named in the config (`Ensembl Gene Id`, `Gene name`, `half_life_PC1`). | yes |
| `gtex_gene_median_tpm.gct.gz` | GTEx Analysis V8, RNA-Seq, "Median gene-level TPM by tissue": `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`, from <https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression>. Save it under the name in the config, or change `paths.gtex`. Used for `MedianExpression_log2`. | no |
| `shap_feature_importance_rankings.csv` | Snapshot of TrunCat's SHAP rankings that defined the candidate set | yes |
| `variants.example.tsv` | Example input: `variant_id, contig, position, refAllele, altAllele, gene, txnames` | yes |

phyloP is no longer used.

Quick check of the GTEx file (a two-line `#1.2` header, then `Name`, `Description`
and one column per tissue):

```bash
zcat gtex_gene_median_tpm.gct.gz | head -3 | cut -c1-120      # run inside Model/TrunKitten/inputs/
```
