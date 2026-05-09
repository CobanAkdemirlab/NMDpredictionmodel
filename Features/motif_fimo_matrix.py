import os
import re
import subprocess
from pathlib import Path
from multiprocessing import Pool

import numpy as np
import pandas as pd
#to install meme.5.5
1) Type the following commands and then follow the instructions
   printed by the "configure" command.

    $ tar zxf meme_VERSION.tar.gz
    $ cd meme_VERSION 
    $ ./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt
# Compile (takes 5-10 minutes)
make
 
# Install
make install
 
# Now the bin directory should exist
ls ~/meme/bin/fimo
 
# Add to PATH
echo 'export PATH=$HOME/meme/bin:$PATH' >> ~/.zshrc
source ~/.zshrc
 
# Test
fimo --version
# CONFIG

WORKDIR = Path("/path/to/location")

SEQUENCE_CSV = WORKDIR / "TOPMed_with_sequences.csv"

FASTA_DIR = WORKDIR / "TOPMed_March_fimo_regions"

MOTIFS = WORKDIR / "memefile.txt"

ATTRACT = WORKDIR / "ATtRACT_db_human.txt"

OUTPUT_DIR = WORKDIR / "TOPMed_March_motif_results"

FIMO_THRESHOLD = 1e-4
FIMO_TIMEOUT = 7200

N_CORES = 6

REGIONS = [
    "utr3_all",
    "utr3_200",
    "ptc_to_ejc",
    "newutr_all",
    "newutr_200",
    "ptc_pm100",
    "ejc_pm100"
]

# MEME PATH

try:
    meme_bin = str(Path.home() / "meme" / "bin")
    meme_lib = str(Path.home() / "meme" / "libexec" / "meme-5.5.5")

    for d in [meme_bin, meme_lib]:
        if os.path.isdir(d) and d not in os.environ.get("PATH", ""):
            os.environ["PATH"] = f"{d}:{os.environ['PATH']}"

except:
    pass
# CHECK FIMO

def check_fimo():

    out = subprocess.run(
        ["fimo", "--version"],
        capture_output=True,
        text=True
    )

    if out.returncode != 0:
        raise RuntimeError("FIMO not found on PATH")

    print("FIMO version:", out.stdout.strip())


# RUN FIMO FOR ONE BATCH
def run_fimo_batch(args):

    region, fasta_path, outdir = args

    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "fimo",
        "--thresh",
        str(FIMO_THRESHOLD),
        "--oc",
        str(outdir),
        str(MOTIFS),
        str(fasta_path)
    ]

    try:

        subprocess.run(
            cmd,
            check=True,
            timeout=FIMO_TIMEOUT,
            capture_output=True,
            text=True
        )

        fimo_file = outdir / "fimo.tsv"

        if fimo_file.exists():
            df = pd.read_csv(fimo_file, sep="\t", comment="#")
            df["region"] = region
            return df

        else:
            return pd.DataFrame()

    except subprocess.TimeoutExpired:

        print("Timeout:", fasta_path)
        return pd.DataFrame()

    except subprocess.CalledProcessError as e:

        print("FIMO error:", e.stderr)
        return pd.DataFrame()
 RUN FIMO FOR ALL REGIONS

def run_fimo():

    print("\nRUNNING FIMO\n")

    tasks = []

    results_dir = OUTPUT_DIR / "fimo_hits"

    for region in REGIONS:

        region_dir = FASTA_DIR / region

        batches = sorted(region_dir.glob("batch_*/seqs.fasta"))

        print(region, "batches:", len(batches))

        for fasta in batches:

            batch_name = fasta.parent.name

            outdir = results_dir / region / batch_name

            tasks.append((region, fasta, outdir))

    with Pool(N_CORES) as pool:

        dfs = pool.map(run_fimo_batch, tasks)

    dfs = [d for d in dfs if not d.empty]

    all_hits = pd.concat(dfs, ignore_index=True)

    hits_file = OUTPUT_DIR / "all_fimo_hits.csv"

    all_hits.to_csv(hits_file, index=False)

    print("Saved motif hits:", hits_file)

    return all_hits

# BUILD MOTIF MATRIX

def build_matrix(hits, n_variants):

    print("\nBUILDING MOTIF MATRIX\n")

    attract = pd.read_csv(ATTRACT, sep="\t")

    id_to_gene = dict(zip(attract["Matrix_id"], attract["Gene_name"]))

    hits["variant_idx"] = (
        hits["sequence_name"]
        .str.replace("seq_", "", regex=False)
        .astype(int)
    )

    def motif_to_gene(m):

        if m in id_to_gene:
            return id_to_gene[m]

        pre = m.split(":")[0]

        if pre in id_to_gene:
            return id_to_gene[pre]

        return ""

    hits["gene"] = hits["motif_id"].apply(motif_to_gene)

    genes = sorted(g for g in hits["gene"].unique() if g)

    matrix = pd.DataFrame(
        0,
        index=np.arange(n_variants),
        columns=[
            f"{r}.{g}"
            for r in REGIONS
            for g in genes
        ],
        dtype=np.int8
    )

    for _, row in hits.iterrows():

        i = row["variant_idx"]

        g = row["gene"]

        r = row["region"]

        if g:

            matrix.loc[i, f"{r}.{g}"] = 1

    matrix_file = OUTPUT_DIR / "topmed_March_motif_presence.csv"

    matrix.to_csv(matrix_file)

    print("Saved motif matrix:", matrix_file)

    return matrix
# MAIN

def main():

    OUTPUT_DIR.mkdir(exist_ok=True)

    check_fimo()

    df = pd.read_csv(SEQUENCE_CSV)

    print("Variants:", len(df))

    hits = run_fimo()

    matrix = build_matrix(hits, len(df))

    print("\nDone.")


if __name__ == "__main__":
    main()
    
