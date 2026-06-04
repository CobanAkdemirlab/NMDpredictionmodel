#!/bin/bash
# ============================================================================
# Title: TOPMed ASE Extraction and Genotype Pipeline
# ============================================================================
# Purpose:
#   1. Extract stop_gained and frameshift variants from TOPMed Freeze 9b
#   2. Remove duplicate/overlapping variants
#   3. Run GATK ASEReadCounter on RNA-seq BAM files
#   4. Extract genotypes from TOPMed BCF for each sample
#   5. Merge ASE positions with sample genotypes
#
# Prerequisites:
#   - bcftools, tabix, GATK installed and in PATH
#   - TOPMed Freeze 9b BCF files (freeze.9b.chrN.pass_and_fail.gtonly.minDP0.bcf)
#   - Sample metadata file with TOR ID → NWD ID mapping
#   - RNA-seq BAM files (*.rna.bam)
#   - Reference genome (GRCh38.fa)
#
# Usage:
#   bash TOPMed_pipeline.sh [config_file]
#
# Output:
#   - Deduplicated PTV VCF: TOPMed_PTVs.dedup.vcf.gz
#   - Per-sample ASE counts: ASE_results/TOR_ID_ASE_counts.table
#   - Per-sample merged genotypes: ASE_genotype/merged_TOR_ID.vcf
# ============================================================================

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

# Load config file if provided, otherwise use defaults
if [ -f "${1:-}" ]; then
    source "$1"
else
    # Default configuration
    INPUT_VCF="/path/to/TOPMed_Freeze9b.vcf.gz"
    METADATA_FILE="/path/to/metadata.txt"
    BAM_DIR="/path/to/bam_files"
    BCF_DIR="/path/to/bcf_files"
    REF_GENOME="/path/to/GRCh38.fa"
    
    # Output directories
    WORK_DIR="$(pwd)"
    OUTPUT_DIR="${WORK_DIR}/TOPMed_output"
    ASE_DIR="${WORK_DIR}/ASE_results"
    GENOTYPE_DIR="${WORK_DIR}/ASE_genotype"
    
    # GATK parameters
    MIN_DEPTH=1
    MIN_MAPPING_QUALITY=255
    MIN_BASE_QUALITY=10
    
    # Parallel processing
    NUM_PARALLEL_JOBS=23
fi

# Create output directories
mkdir -p "${OUTPUT_DIR}" "${ASE_DIR}" "${GENOTYPE_DIR}"

echo "========================================"
echo "TOPMed ASE Extraction and Genotype Pipeline"
echo "========================================"
echo "Work directory: ${WORK_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# ============================================================================
# PART 1: Extract and Deduplicate Variants
# ============================================================================

echo "[STEP 1] Extracting stop_gained and frameshift variants..."

bcftools view \
    --drop-genotypes \
    -i 'INFO/ANN ~ "frameshift_variant" || INFO/ANN ~ "stop_gained"' \
    "${INPUT_VCF}" \
    -Oz \
    -o "${OUTPUT_DIR}/TOPMed_PTVs.vcf.gz"

echo "✓ Extracted PTV VCF: ${OUTPUT_DIR}/TOPMed_PTVs.vcf.gz"

echo ""
echo "[STEP 2] Removing duplicate/overlapping variants..."

bcftools norm \
    --rm-dup all \
    "${OUTPUT_DIR}/TOPMed_PTVs.vcf.gz" \
    -Oz \
    -o "${OUTPUT_DIR}/TOPMed_PTVs.dedup.vcf.gz"

echo "✓ Deduplicated PTV VCF: ${OUTPUT_DIR}/TOPMed_PTVs.dedup.vcf.gz"

echo ""
echo "[STEP 2b] Indexing deduplicated VCF..."

tabix -p vcf "${OUTPUT_DIR}/TOPMed_PTVs.dedup.vcf.gz"

echo "✓ VCF indexed"

# ============================================================================
# PART 2: Extract ASE Counts from RNA-seq BAM Files
# ============================================================================

echo ""
echo "[STEP 3] Extracting ASE counts using GATK ASEReadCounter..."

# Find all RNA-seq BAM files
bam_files=($(find "${BAM_DIR}" -name "*.rna.bam" -type f))

if [ ${#bam_files[@]} -eq 0 ]; then
    echo "Error: No BAM files found in ${BAM_DIR}"
    exit 1
fi

echo "Found ${#bam_files[@]} RNA-seq BAM files"

for bam_file in "${bam_files[@]}"; do
    sample_name=$(basename "$bam_file" .rna.bam)
    output_table="${ASE_DIR}/${sample_name}_ASE_counts.table"
    
    echo "  Processing: ${sample_name}"
    
    gatk ASEReadCounter \
        -R "${REF_GENOME}" \
        -I "$bam_file" \
        -V "${OUTPUT_DIR}/TOPMed_PTVs.dedup.vcf.gz" \
        -O "$output_table" \
        --min-depth ${MIN_DEPTH} \
        --min-mapping-quality ${MIN_MAPPING_QUALITY} \
        --min-base-quality ${MIN_BASE_QUALITY}
    
    echo "    ✓ ${output_table}"
done

# ============================================================================
# PART 3: Extract Genotypes for Each Sample
# ============================================================================

echo ""
echo "[STEP 4] Extracting genotypes from TOPMed BCF for each sample..."

# Function to process a single sample
process_sample() {
    local sample_tor="$1"
    local metadata_file="$2"
    local bcf_dir="$3"
    local genotype_dir="$4"
    local ase_dir="$5"
    local col_num="$6"
    local num_jobs="$7"
    
    # Extract NWD ID from metadata
    local nwd_id=$(grep "${sample_tor}" "${metadata_file}" | awk '{print $6}')
    
    if [ -z "$nwd_id" ]; then
        echo "Error: Could not find ${sample_tor} in metadata"
        return 1
    fi
    
    echo "  Processing ${sample_tor} (NWD: ${nwd_id})"
    
    # Create temporary working directory
    local work_dir="/tmp/${sample_tor}_processing"
    mkdir -p "${work_dir}"
    cd "${work_dir}"
    
    # Process each chromosome
    process_chr() {
        local chr="$1"
        local ase_csv="${ase_dir}/${sample_tor}_ASE_counts.table"
        
        if [ ! -f "$ase_csv" ]; then
            echo "Warning: ASE file not found: $ase_csv"
            return
        fi
        
        # Extract positions from ASE CSV for this chromosome
        awk -F'\t' -v target_chr="${chr}" 'NR > 1 && $1 == target_chr {print $1, $2}' OFS='\t' \
            "${ase_csv}" > "${chr}_positions.txt"
        
        if [ ! -s "${chr}_positions.txt" ]; then
            return
        fi
        
        # Extract genotypes from BCF
        bcftools view \
            -R "${chr}_positions.txt" \
            "${bcf_dir}/freeze.9b.${chr}.pass_and_fail.gtonly.minDP0.bcf" \
            -Oz \
            -o "${chr}_genotypes.vcf.gz"
        
        # Extract columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE_GENOTYPE
        bcftools view "${chr}_genotypes.vcf.gz" | \
            grep -v '##' | \
            cut -f "1,2,3,4,5,6,7,9,${col_num}" > "${chr}_genotype.txt"
        
        # Extract ASE positions
        awk -F'\t' -v target_chr="${chr}" 'NR == 1 || ($1 == target_chr) {print}' \
            "${ase_csv}" > "${chr}_positions.tsv"
        
        # Merge genotypes with ASE data
        awk 'BEGIN {FS=OFS="\t"} 
            FNR==NR {
                genotype[$1, $2] = $9
                next
            } 
            FNR == 1 {
                print $0, "'"${nwd_id}"'"
                next
            } 
            ($1, $2) in genotype {
                print $0, genotype[$1, $2]
            }' \
            "${chr}_genotype.txt" "${chr}_positions.tsv" > "merged_${chr}_positions.vcf"
    }
    
    export -f process_chr
    
    # Run chromosome processing in parallel
    parallel -j "${num_jobs}" process_chr ::: \
        $(seq 1 22 | sed 's/^/chr/') chrX chrY chrM 2>/dev/null || true
    
    # Merge all chromosomes
    if ls merged_chr*_positions.vcf 1>/dev/null 2>&1; then
        merged_files=$(ls merged_chr*_positions.vcf | sort -V)
        
        (head -n 1 $(echo "$merged_files" | head -n 1) && \
         tail -n +2 -q $(echo "$merged_files")) > "merged_${sample_tor}.vcf"
        
        # Copy to final output directory
        mv "merged_${sample_tor}.vcf" "${genotype_dir}/"
        echo "    ✓ ${genotype_dir}/merged_${sample_tor}.vcf"
    fi
    
    # Cleanup
    cd /
    rm -rf "${work_dir}"
}

export -f process_sample

# Get column number for each sample and process
while IFS= read -r metadata_line; do
    sample_tor=$(echo "$metadata_line" | awk '{print $3}')
    nwd_id=$(echo "$metadata_line" | awk '{print $6}')
    
    if [ -z "$sample_tor" ] || [ -z "$nwd_id" ]; then
        continue
    fi
    
    # Find column number for this sample in BCF
    echo "  Finding column for ${sample_tor}..."
    col_num=$(bcftools view -h "${BCF_DIR}/freeze.9b.chr1.pass_and_fail.gtonly.minDP0.bcf" | \
        awk -F'\t' -v nwd="${nwd_id}" '{for(i=1; i<=NF; i++) if($i == nwd) {print i; exit}}')
    
    if [ -z "$col_num" ]; then
        echo "    Error: Could not find column for ${nwd_id}"
        continue
    fi
    
    echo "    Column: ${col_num}"
    
    process_sample "${sample_tor}" "${METADATA_FILE}" "${BCF_DIR}" "${GENOTYPE_DIR}" "${ASE_DIR}" "${col_num}" "${NUM_PARALLEL_JOBS}"
    
done < <(tail -n +2 "${METADATA_FILE}")

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "========================================"
echo "Pipeline Complete!"
echo "========================================"
echo ""
echo "Output files:"
echo "  • Deduplicated VCF: ${OUTPUT_DIR}/TOPMed_PTVs.dedup.vcf.gz"
echo "  • ASE counts: ${ASE_DIR}/*.table"
echo "  • Merged genotypes: ${GENOTYPE_DIR}/merged_*.vcf"
echo ""
echo "Next step: Run TOPMed_merge_and_annotate.R"
echo ""
