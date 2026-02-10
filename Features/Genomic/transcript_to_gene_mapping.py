import pandas as pd
import numpy as np
import gzip
import re

# --- Paths ---
csv_path = "~/TOPMed_df_snv_feb2026.csv"
gtf_path = "~/gencode.v26.primary_assembly.annotation.gtf.gz"

# --- Step 1: Load your dataset and get unique transcript names ---
df = pd.read_csv(csv_path)
txnames = df["txnames"].dropna().unique().tolist()

print(f"Loaded {len(txnames):,} unique transcript IDs from txnames column")

# --- Step 2: Parse GTF to extract transcript_id ↔ gene_id mapping ---
mapping = []

with gzip.open(gtf_path, "rt") as gtf:
    for line in gtf:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue
        attr = fields[8]

        # Extract transcript_id and gene_id
        transcript_match = re.search(r'transcript_id "([^"]+)"', attr)
        gene_match = re.search(r'gene_id "([^"]+)"', attr)

        if transcript_match and gene_match:
            transcript_id = transcript_match.group(1)
            gene_id = gene_match.group(1)
            mapping.append((transcript_id, gene_id))

# Convert to DataFrame and drop duplicates
tx2gene = pd.DataFrame(mapping, columns=["transcript_id", "gene_id"]).drop_duplicates()

print(f"Extracted {len(tx2gene):,} transcript–gene pairs from GTF")

# --- Step 3: Filter to only transcripts in your dataset ---
# Normalize both transcript ID lists by removing version suffixes
tx2gene["base_tx_id"] = tx2gene["transcript_id"].str.replace(r"\.\d+$", "", regex=True)
df["base_tx_id"] = df["txnames"].astype(str).str.replace(r"\.\d+$", "", regex=True)

# Re-match using normalized IDs
subset_map = tx2gene[tx2gene["base_tx_id"].isin(df["base_tx_id"])]
print(f"Matched {len(subset_map):,} of {len(df['base_tx_id'].unique()):,} transcripts after normalization")

# Save normalized mapping
subset_map[["transcript_id", "gene_id"]].drop_duplicates().to_csv(
    "TOPMed_tx_to_gene_mapping_from_gencode_v26_normalized.csv", index=False
)
print("✅ Saved normalized mapping: TOPMed_tx_to_gene_mapping_from_gencode_v26_normalized.csv")
