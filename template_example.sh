#!/usr/bin/env bash

# ==============================================================================
# Script Name :  run_qc_pipeline.sh
# Description :  Runs FastQC and fastp on paired-end genomic data.
# Author      :  Songphon Sutthitthasakul
# Updated Date:  25-01-2026
# ==============================================================================

set -euo pipefail

# 1. VARIABLE SETTING
# Use readonly for variables that should not change
readonly THREADS=8
readonly MEMORY=20 # Gb unit
readonly IN_DIR="./data/raw"
readonly OUT_DIR="./data/processed"
readonly REPORT_DIR="./reports"

# 2. FUNCTIONS
# --- Function: Run FastQC ---
run_fastqc() {
  local input_file=$1
  local output_path=$2

  echo "LOG: Running FastQC on ${input_file}..."
  mkdir -p "${output_path}"

  fastqc -t "${THREADS}" "${input_file}" -o "${output_path}"
}

# --- Function: Run fastp (Trimming & Filtering) ---
run_fastp() {
  local r1_in=$1
  local r2_in=$2
  local sample_id=$3

  echo "LOG: Running fastp for sample ${sample_id}..."
  mkdir -p "${OUT_DIR}" "${REPORT_DIR}"

  fastp \
    --thread "${THREADS}" \
    --in1 "${r1_in}" \
    --in2 "${r2_in}" \
    --out1 "${OUT_DIR}/${sample_id}_R1_clean.fastq.gz" \
    --out2 "${OUT_DIR}/${sample_id}_R2_clean.fastq.gz" \
    --json "${REPORT_DIR}/${sample_id}.json" \
    --html "${REPORT_DIR}/${sample_id}.html" \
    --detect_adapter_for_pe \
    --trim_poly_g
}

# 3. RUNNING (Main Logic)
main() {
  # Ensure directories exist
  mkdir -p "${OUT_DIR}" "${REPORT_DIR}"

  # Loop through R1 files and find matching R2
  for r1_file in "${IN_DIR}"/*_R1.fastq.gz; do
    # Extract Sample ID (e.g., sample01_R1.fastq.gz -> sample01)
    local base_name
    base_name=$(basename "${r1_file}" _R1.fastq.gz)
    local r2_file="${IN_DIR}/${base_name}_R2.fastq.gz"

    # 1. Initial QC
    run_fastqc "${r1_file}" "${REPORT_DIR}/fastqc_raw"
    run_fastqc "${r2_file}" "${REPORT_DIR}/fastqc_raw"

    # 2. Processing with fastp
    run_fastp "${r1_file}" "${r2_file}" "${base_name}"

    # 3. Post-processing QC
    run_fastqc "${OUT_DIR}/${base_name}_R1_clean.fastq.gz" "${REPORT_DIR}/fastqc_clean"
  done

  echo "SUCCESS: Pipeline completed."
}

# Execution
main "$@"
