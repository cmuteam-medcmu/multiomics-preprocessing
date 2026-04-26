#!/bin/bash
set -euo pipefail

# =============================================================================
# Script Name :  summary.sh
# Description :  Run calc_region_depth + MultiQC + QC summary report
# Usage       :  bash summary.sh <OUT_DIR> <QC_SCRIPT> [SAMPLE_REGEX]
# =============================================================================

# --- Arguments ---
OUT_DIR="${1:?Usage: bash summary.sh <OUT_DIR> <QC_SCRIPT> [SAMPLE_REGEX]}"
QC_SCRIPT="${2:?Usage: bash summary.sh <OUT_DIR> <QC_SCRIPT> [SAMPLE_REGEX]}"
DEFAULT_REGEX='^([^_]+(?:_[^_]+){5})'
SAMPLE_REGEX="${3:-$DEFAULT_REGEX}"

REPORT_DIR="${OUT_DIR}/reports"
PREPROCESS_DIR="${OUT_DIR}/02_PREPROCESS"
DEPTH_TSV="${REPORT_DIR}/mosdepth_region_depth.tsv"

# --- 1. Calculate chr1-22 region depth for each sample ---
echo "LOG: Calculating chr1-22 region depth for all samples..."

# Create header
echo -e "Sample\tMean_Region_Depth_chr1_22\tMedian_Region_Depth_chr1_22" > "${DEPTH_TSV}"

# Loop over all regions.bed.gz files in preprocess subdirectories
for regions_file in "${PREPROCESS_DIR}"/*/*.regions.bed.gz; do
    [ -f "${regions_file}" ] || continue

    sample_id=$(basename "${regions_file}" .regions.bed.gz)
    echo "  Processing ${sample_id}..."

    # Extract chr1-22 depth values, sort numerically, then compute mean+median
    stats=$(zcat "${regions_file}" \
      | awk '$1 ~ /^chr[0-9]+$/ { n=substr($1,4)+0; if(n>=1 && n<=22) print $4 }' \
      | sort -g \
      | awk '{sum+=$1; vals[NR]=$1}
        END {
          if(NR==0){printf "NA\tNA\n"; exit}
          mean=sum/NR
          if(NR%2==1) median=vals[int(NR/2)+1]
          else median=(vals[NR/2]+vals[NR/2+1])/2
          printf "%.4f\t%.4f\n", mean, median
        }')

    echo -e "${sample_id}\t${stats}" >> "${DEPTH_TSV}"
done

echo "LOG: Region depth TSV written to ${DEPTH_TSV}"

# --- 2. MultiQC ---
ml purge
ml multiqc

echo "LOG: Running MultiQC on ${OUT_DIR}..."
multiqc "${OUT_DIR}" -o "${OUT_DIR}" --force

# --- 3. QC Summary Report ---
echo "LOG: Running QC summary script..."
python "${QC_SCRIPT}" \
    -i "${OUT_DIR}/multiqc_data/multiqc_data.json" \
    -o "${OUT_DIR}/QC_SUMMARY.xlsx" \
    -p "${SAMPLE_REGEX}" \
    -d "${DEPTH_TSV}"

echo "SUCCESS: Summary complete. Output: ${OUT_DIR}/QC_SUMMARY.xlsx"
