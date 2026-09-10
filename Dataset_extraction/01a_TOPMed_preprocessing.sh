#!/usr/bin/env bash

# ==============================================================================
# TOPMed ASE extraction and genotype-matching pipeline
#
# Purpose:
#   1. Extract candidate stop-gained and frameshift variants from TOPMed
#      Freeze 9b
#   2. Remove duplicate variant records
#   3. Run GATK ASEReadCounter on matched RNA-seq BAM files
#   4. Retrieve donor genotypes at ASE-observed positions from TOPMed Freeze 9b
#      genotype BCFs
#   5. Merge ASE counts with donor genotypes
#
# Input:
#   - TOPMed Freeze 9b annotated VCF
#   - TOPMed Freeze 9b per-chromosome genotype BCFs
#   - Matched RNA-seq BAM files
#   - Sample metadata mapping TOR IDs to NWD IDs
#   - GRCh38 reference FASTA
#
# Output:
#   - TOPMed_PTVs.dedup.vcf.gz
#   - Per-sample ASEReadCounter tables
#   - Per-sample merged ASE/genotype files: merged_TOR*.vcf
#
# Next step:
#   01b_TOPMed_extraction.R
# ==============================================================================

set -euo pipefail


# ==============================================================================
# Configuration
# ==============================================================================

if [ -f "${1:-}" ]; then
    source "$1"
else
    INPUT_VCF="/path/to/TOPMed_Freeze9b.vcf.gz"
    METADATA_FILE="/path/to/metadata.tsv"
    BAM_DIR="/path/to/bam_files"
    BCF_DIR="/path/to/bcf_files"
    REF_GENOME="/path/to/GRCh38.fa"

    WORK_DIR="$(pwd)"
    OUTPUT_DIR="${WORK_DIR}/TOPMed_output"
    ASE_DIR="${WORK_DIR}/ASE_results"
    GENOTYPE_DIR="${WORK_DIR}/ASE_genotype"

    MIN_DEPTH=1
    MIN_MAPPING_QUALITY=255
    MIN_BASE_QUALITY=10
fi

mkdir -p \
    "${OUTPUT_DIR}" \
    "${ASE_DIR}" \
    "${GENOTYPE_DIR}"


# ==============================================================================
# STEP 1: Extract candidate PTV sites
# ==============================================================================

echo "[STEP 1] Extracting candidate stop-gained and frameshift variants..."

bcftools view \
    --drop-genotypes \
    -i 'INFO/ANN ~ "frameshift_variant" || INFO/ANN ~ "stop_gained"' \
    "${INPUT_VCF}" \
    -Oz \
    -o "${OUTPUT_DIR}/TOPMed_PTVs.vcf.gz"


# ==============================================================================
# STEP 2: Remove duplicate variant records
# ==============================================================================

echo "[STEP 2] Removing duplicate variant records..."

bcftools norm \
    --rm-dup all \
    "${OUTPUT_DIR}/TOPMed_PTVs.vcf.gz" \
    -Oz \
    -o "${OUTPUT_DIR}/TOPMed_PTVs.dedup.vcf.gz"

tabix -f -p vcf \
    "${OUTPUT_DIR}/TOPMed_PTVs.dedup.vcf.gz"


# ==============================================================================
# STEP 3: ASE extraction
# ==============================================================================

echo "[STEP 3] Running GATK ASEReadCounter..."

find "${BAM_DIR}" \
    -type f \
    -name "*.rna.bam" \
    -print0 |
while IFS= read -r -d '' bam_file; do

    sample_tor=$(basename "${bam_file}" .rna.bam)

    echo "  Processing ${sample_tor}"

    gatk ASEReadCounter \
        -R "${REF_GENOME}" \
        -I "${bam_file}" \
        -V "${OUTPUT_DIR}/TOPMed_PTVs.dedup.vcf.gz" \
        -O "${ASE_DIR}/${sample_tor}_ASE_counts.table" \
        --min-depth "${MIN_DEPTH}" \
        --min-mapping-quality "${MIN_MAPPING_QUALITY}" \
        --min-base-quality "${MIN_BASE_QUALITY}"

done


# ==============================================================================
# STEP 4: Retrieve donor genotypes and merge with ASE
# ==============================================================================

echo "[STEP 4] Matching ASE observations with donor genotypes..."

while IFS=$'\t' read -r sample_tor nwd_id; do

    # Skip header if present
    if [ "${sample_tor}" = "TOR_ID" ]; then
        continue
    fi

    ase_file="${ASE_DIR}/${sample_tor}_ASE_counts.table"

    if [ ! -f "${ase_file}" ]; then
        echo "  Warning: ASE file not found for ${sample_tor}"
        continue
    fi

    echo "  Processing ${sample_tor} (${nwd_id})"

    sample_tmp=$(mktemp -d)

    for chr in \
        chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
        chr21 chr22 chrX chrY; do

        positions_file="${sample_tmp}/${chr}_positions.txt"

        awk \
            -F'\t' \
            -v target_chr="${chr}" \
            'NR > 1 && $1 == target_chr {print $1 "\t" $2}' \
            "${ase_file}" \
            > "${positions_file}"

        if [ ! -s "${positions_file}" ]; then
            continue
        fi

        bcf_file="${BCF_DIR}/freeze.9b.${chr}.pass_and_fail.gtonly.minDP0.bcf"

        if [ ! -f "${bcf_file}" ]; then
            echo "    Warning: missing ${bcf_file}"
            continue
        fi

        # Extract this donor's genotype at ASE-observed positions
        bcftools query \
            -s "${nwd_id}" \
            -R "${positions_file}" \
            -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
            "${bcf_file}" \
            > "${sample_tmp}/${chr}_genotypes.tsv"

    done

    # Combine chromosome-level genotype files
    cat "${sample_tmp}"/chr*_genotypes.tsv 2>/dev/null \
        > "${sample_tmp}/genotypes_all.tsv" || true

    # Merge genotype with ASE by chromosome + position
    awk '
        BEGIN {FS=OFS="\t"}

        FNR==NR {
            genotype[$1 SUBSEP $2] = $5
            next
        }

        FNR==1 {
            print $0, "GT"
            next
        }

        ($1 SUBSEP $2) in genotype {
            print $0, genotype[$1 SUBSEP $2]
        }
    ' \
        "${sample_tmp}/genotypes_all.tsv" \
        "${ase_file}" \
        > "${GENOTYPE_DIR}/merged_${sample_tor}.vcf"

    rm -rf "${sample_tmp}"

done < "${METADATA_FILE}"


# ==============================================================================
# Summary
# ==============================================================================

echo ""
echo "TOPMed preprocessing complete."
echo ""
echo "Candidate PTV VCF:"
echo "  ${OUTPUT_DIR}/TOPMed_PTVs.dedup.vcf.gz"
echo ""
echo "ASEReadCounter output:"
echo "  ${ASE_DIR}/*_ASE_counts.table"
echo ""
echo "Merged ASE/genotype files:"
echo "  ${GENOTYPE_DIR}/merged_TOR*.vcf"
echo ""
echo "Next step:"
echo "  Rscript 01b_TOPMed_extraction.R"
