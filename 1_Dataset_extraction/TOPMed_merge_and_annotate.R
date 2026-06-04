# ============================================================================
# Title: TOPMed Heterozygous PTC Variant Merging and ANNOVAR Annotation
# ============================================================================
# Purpose:
#   1. Merge all per-sample genotype VCF files (merged_TOR*.vcf)
#   2. Filter for heterozygous (0/1) protein-truncating variants
#   3. Prepare merged data for ANNOVAR annotation
#   4. Run ANNOVAR functional annotation
#   5. Generate summary statistics
#
# Prerequisites:
#   - Run TOPMed_pipeline.sh to generate merged_TOR*.vcf files
#   - ANNOVAR installed and configured
#
# Input:
#   - Per-sample merged VCF files: ASE_genotype/merged_TOR*.vcf
#
# Output:
#   - Merged heterozygous variants: NMD_TOPMed/merged_het_variants.txt
#   - VCF for ANNOVAR: NMD_TOPMed/TOPMed_stopgain_frameshift.vcf
#   - ANNOVAR input: NMD_TOPMed/TOPMed_stopgain_frameshift.avinput
#   - ANNOVAR annotation: NMD_TOPMed/TOPMed_stopgain_frameshift_annotated.*
# ============================================================================

library(data.table)
library(dplyr)

# ============================================================================
# Configuration
# ============================================================================

CONFIG <- list(
  input_dir = "~/ASE_genotype",
  output_dir = "~/NMD_TOPMed",
  annovar_path = "~/annovar",
  annovar_buildver = "hg38",
  annovar_db = "ensGene",
  vcf_pattern = "^merged_.*\\.vcf$"
)

# ============================================================================
# Create output directory
# ============================================================================

dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# PART 1: Find and Read All Per-Sample VCF Files
# ============================================================================

message("=" %*% 70)
message("TOPMed Heterozygous PTC Variant Merging and ANNOVAR Annotation")
message("=" %*% 70)

setwd(CONFIG$input_dir)

# List all merged sample VCF files
vcf_files <- list.files(
  ".",
  pattern = CONFIG$vcf_pattern,
  full.names = FALSE
)

if (length(vcf_files) == 0) {
  stop("No merged VCF files found in: ", CONFIG$input_dir)
}

message("\n[STEP 1] Reading per-sample VCF files")
message("Found ", length(vcf_files), " sample files:")
for (f in vcf_files) message("  - ", f)

# ============================================================================
# PART 2: Read and Filter VCF Files
# ============================================================================

message("\n[STEP 2] Filtering for heterozygous (0/1) variants...")

het_variants_list <- lapply(seq_along(vcf_files), function(idx) {
  file <- vcf_files[idx]
  sample_name <- sub("^merged_|.vcf$", "", file)
  
  message("  [", idx, "/", length(vcf_files), "] ", sample_name, " ", appendLF = FALSE)
  
  tryCatch({
    # Read VCF file
    df <- data.table::fread(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    if (nrow(df) == 0) {
      message("(empty)")
      return(NULL)
    }
    
    # Get the sample column (should be the last column with genotypes)
    last_col_idx <- ncol(df)
    sample_col_name <- colnames(df)[last_col_idx]
    
    # Filter for heterozygous genotypes (0/1)
    het_df <- df[df[[last_col_idx]] == "0/1", ]
    
    if (nrow(het_df) == 0) {
      message("(no het variants)")
      return(NULL)
    }
    
    # Add sample identifier
    het_df[, sample_id := sample_name]
    
    message("(", nrow(het_df), " variants)")
    return(het_df)
    
  }, error = function(e) {
    message("(ERROR: ", conditionMessage(e), ")")
    return(NULL)
  })
})

# Remove NULL entries
het_variants_list <- Filter(Negate(is.null), het_variants_list)

if (length(het_variants_list) == 0) {
  stop("No heterozygous variants found in any sample VCF files")
}

# Combine all samples into single data.table
het_get <- data.table::rbindlist(het_variants_list, fill = TRUE)

message("\n✓ Total heterozygous variants: ", nrow(het_get))
message("✓ Unique samples: ", n_distinct(het_get$sample_id))

# ============================================================================
# PART 3: Save Merged Heterozygous Variants
# ============================================================================

message("\n[STEP 3] Saving merged heterozygous variants...")

merged_output <- file.path(CONFIG$output_dir, "merged_het_variants.txt")

data.table::fwrite(
  het_get,
  file = merged_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message("✓ Saved to: ", merged_output)

# ============================================================================
# PART 4: Standardize and Prepare for VCF Output
# ============================================================================

message("\n[STEP 4] Preparing data for ANNOVAR annotation...")

# Ensure VCF standard columns exist
vcf_required_cols <- c("CHROM", "POS", "ID", "REF", "ALT")

for (col in c("CHROM", "POS")) {
  if (!(col %in% colnames(het_get))) {
    # Try alternative column names from VCF standard
    alt_name <- tolower(col)
    if (alt_name %in% tolower(colnames(het_get))) {
      actual_col <- colnames(het_get)[grep(paste0("^", alt_name, "$"), colnames(het_get), ignore.case = TRUE)][1]
      if (!is.na(actual_col)) {
        het_get[[col]] <- het_get[[actual_col]]
      }
    }
  }
}

# Ensure ID column exists
if (!("ID" %in% colnames(het_get))) {
  het_get[, ID := "."]
}

# Get REF and ALT from standard VCF columns
if (!("REF" %in% colnames(het_get))) {
  ref_col <- colnames(het_get)[grep("^REF|^ref", colnames(het_get))][1]
  if (!is.na(ref_col)) {
    het_get[, REF := get(ref_col)]
  } else {
    message("Warning: REF column not found")
  }
}

if (!("ALT" %in% colnames(het_get))) {
  alt_col <- colnames(het_get)[grep("^ALT|^alt", colnames(het_get))][1]
  if (!is.na(alt_col)) {
    het_get[, ALT := get(alt_col)]
  } else {
    message("Warning: ALT column not found")
  }
}

# Add VCF format columns if missing
if (!("QUAL" %in% colnames(het_get))) het_get[, QUAL := "."]
if (!("FILTER" %in% colnames(het_get))) het_get[, FILTER := "PASS"]
if (!("INFO" %in% colnames(het_get))) het_get[, INFO := "."]
if (!("FORMAT" %in% colnames(het_get))) het_get[, FORMAT := "GT"]

# Select columns in VCF order
vcf_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

# Check which columns exist
vcf_cols_present <- intersect(vcf_cols, colnames(het_get))

if (length(vcf_cols_present) < 5) {
  stop("Missing required VCF columns. Found: ", paste(vcf_cols_present, collapse = ", "))
}

vcf_output <- het_get[, ..vcf_cols_present]

message("✓ Data prepared for ANNOVAR")
message("  Columns: ", paste(colnames(vcf_output), collapse = ", "))

# ============================================================================
# PART 5: Write VCF File for ANNOVAR
# ============================================================================

message("\n[STEP 5] Writing VCF file for ANNOVAR...")

vcf_path <- file.path(CONFIG$output_dir, "TOPMed_stopgain_frameshift.vcf")

data.table::fwrite(
  vcf_output,
  file = vcf_path,
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)

message("✓ VCF file: ", vcf_path)
message("  Total variants: ", nrow(vcf_output))

# ============================================================================
# PART 6: Run ANNOVAR Annotation
# ============================================================================

message("\n[STEP 6] Running ANNOVAR annotation...")

# Step 6a: Convert VCF to ANNOVAR format
message("\n  [6a] Converting VCF to ANNOVAR format...")

avinput_path <- file.path(CONFIG$output_dir, "TOPMed_stopgain_frameshift.avinput")

convert_cmd <- sprintf(
  "perl %s/convert2annovar.pl -format vcf4 %s > %s 2>&1",
  CONFIG$annovar_path,
  vcf_path,
  avinput_path
)

message("  Command: ", convert_cmd)
result <- system(convert_cmd, wait = TRUE)

if (result != 0) {
  warning("ANNOVAR conversion may have failed (exit code: ", result, ")")
} else if (!file.exists(avinput_path) || file.size(avinput_path) == 0) {
  warning("ANNOVAR output file is empty or missing")
} else {
  message("  ✓ ANNOVAR input: ", avinput_path)
}

# Step 6b: Run ANNOVAR functional annotation
message("\n  [6b] Running functional annotation with ANNOVAR...")

annotated_prefix <- file.path(CONFIG$output_dir, "TOPMed_stopgain_frameshift_annotated")

annotate_cmd <- sprintf(
  "perl %s/annotate_variation.pl -buildver %s -out %s -dbtype %s %s %s/humandb/ 2>&1",
  CONFIG$annovar_path,
  CONFIG$annovar_buildver,
  annotated_prefix,
  CONFIG$annovar_db,
  avinput_path,
  CONFIG$annovar_path
)

message("  Command: ", annotate_cmd)
result <- system(annotate_cmd, wait = TRUE)

if (result != 0) {
  warning("ANNOVAR annotation may have failed (exit code: ", result, ")")
  message("  Please check ANNOVAR installation and paths")
} else {
  annotated_files <- c(
    sprintf("%s.%s_multianno.txt", annotated_prefix, CONFIG$annovar_buildver),
    sprintf("%s.%s_multianno.vcf", annotated_prefix, CONFIG$annovar_buildver),
    sprintf("%s.%s_multianno.xlsx", annotated_prefix, CONFIG$annovar_buildver)
  )
  
  existing_files <- annotated_files[file.exists(annotated_files)]
  
  if (length(existing_files) > 0) {
    message("  ✓ ANNOVAR annotation complete")
    for (f in existing_files) {
      message("    - ", basename(f))
    }
  } else {
    message("  Warning: No ANNOVAR output files found")
  }
}

# ============================================================================
# PART 7: Summary Report
# ============================================================================

message("\n" %*% 70)
message("ANALYSIS COMPLETE")
message("=" %*% 70)

message("\n[SUMMARY STATISTICS]")
message("  • Input VCF files processed: ", length(vcf_files))
message("  • Total heterozygous variants: ", nrow(het_get))
message("  • Unique samples: ", n_distinct(het_get$sample_id))
message("  • Variants written to VCF: ", nrow(vcf_output))

message("\n[OUTPUT FILES]")
message("  Directory: ", CONFIG$output_dir)
message("")
message("  1. merged_het_variants.txt")
message("     Complete merged dataset with all variants and metadata")
message("")
message("  2. TOPMed_stopgain_frameshift.vcf")
message("     VCF format file (input for ANNOVAR)")
message("")
message("  3. TOPMed_stopgain_frameshift.avinput")
message("     ANNOVAR format input file")
message("")
message("  4. TOPMed_stopgain_frameshift_annotated.", CONFIG$annovar_buildver, "_multianno.txt")
message("     ANNOVAR functional annotation results (text format)")
message("")
message("  5. TOPMed_stopgain_frameshift_annotated.", CONFIG$annovar_buildver, "_multianno.vcf")
message("     ANNOVAR annotated VCF")
message("")
message("  6. TOPMed_stopgain_frameshift_annotated.", CONFIG$annovar_buildver, "_multianno.xlsx")
message("     ANNOVAR results in Excel format")

message("\n[CONFIGURATION]")
message("  Input directory: ", CONFIG$input_dir)
message("  Output directory: ", CONFIG$output_dir)
message("  ANNOVAR build: ", CONFIG$annovar_buildver)
message("  ANNOVAR database: ", CONFIG$annovar_db)

message("\n" %*% 70)
message("Script completed successfully!")
message("=" %*% 70)
