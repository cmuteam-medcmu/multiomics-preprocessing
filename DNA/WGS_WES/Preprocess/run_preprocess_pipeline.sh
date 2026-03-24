#!/bin/bash

#SBATCH --account=o250022               
#SBATCH --job-name=fq2bam  
#SBATCH --partition=gpu
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=40
#SBATCH --gpus=1
#SBATCH --mem=150G
#SBATCH --output=preprocess_%j.out          
#SBATCH --error=preprocess_%j.err

# ==============================================================================
# Script Name :  run_preprocess_pipeline.sh
# Description :  Runs FastQC, Fastp, and FQ2BAM on paired-end genomic data.
# Author      :  Songphon Sutthitthasakul
# Updated Date:  11-Mar-2026
# ==============================================================================

set -euo pipefail

# 1. VARIABLE SETTING
readonly THREADS=40
readonly MEMORY=200 # Gb unit
readonly OUT_DIR="/project/o250022_cfOSTEO/Cell-line/20260320_Tissue_WES"
readonly SAMPLES_SHEET="/project/o250022_cfOSTEO/Cell-line/20260320_Tissue_WES/00_RAW_FASTQ/sample_sheet.csv"
readonly QC_SCRIPT="/project/o250022_cfOSTEO/Cell-line/script/qc_script.py"
readonly REPORT_DIR="${OUT_DIR}/reports"

# 2. FUNCTIONS
# --- Function: Run FastQC ---
run_fastqc() {
  ml purge
  ml fastqc

  local input_r1_file=$1
  local input_r2_file=$2
  local output_path=$3

  echo "LOG: Running FastQC on ${input_r1_file}..."
  echo "LOG: Running FastQC on ${input_r2_file}..."
  mkdir -p "${output_path}"

  fastqc -t 2 "${input_r1_file}" "${input_r2_file}" -o "${output_path}"
}

# --- Function: Run fastp (Trimming & Filtering) ---
run_fastp() {
  ml purge
  ml fastp

  local r1_in=$1
  local r2_in=$2
  local sample_id=$3
  local output_dir=$4

  echo "LOG: Running Fastp for sample ${sample_id}..."
  mkdir -p "${output_dir}" "${REPORT_DIR}"

  fastp \
    --thread "${THREADS}" \
    --in1 "${r1_in}" \
    --in2 "${r2_in}" \
    --out1 "${output_dir}/${sample_id}_R1_clean.fastq.gz" \
    --out2 "${output_dir}/${sample_id}_R2_clean.fastq.gz" \
    --json "${REPORT_DIR}/${sample_id}.json" \
    --html "${REPORT_DIR}/${sample_id}.html" \
    --detect_adapter_for_pe \
    --trim_poly_g
}

# --- Function: Run FQ2BAM ---
run_fq2bam() {
  ml purge
  ml fastqc/0.12.1
  ml apptainer
  ml samtools

  local sample_id=$1
  local input_dir=$2
  local output_dir="${3}/${sample_id}"
  local FASTA_FILE="/common/db/human_ref/hg38/parabricks/Homo_sapiens_assembly38.fasta"

  echo "LOG: Running FQ2BAM for sample ${sample_id}..."
  mkdir -p "${output_dir}" "${REPORT_DIR}"

  # FQ2BAM
  apptainer exec --nv \
    -B /common/db/human_ref/hg38/parabricks:/Ref_hg38 \
    -B /common/db/human_ref/hg38/v0:/Ref_hg38_v0 \
    -B ${input_dir}:/FASTQdir \
    -B ${output_dir}:/Workdir \
    /common/sif/clara-parabricks/4.4.0.sif pbrun \
      fq2bam \
      --ref /Ref_hg38/Homo_sapiens_assembly38.fasta \
      --in-fq /FASTQdir/${sample_id}_R1_clean.fastq.gz /FASTQdir/${sample_id}_R2_clean.fastq.gz \
      --knownSites /Ref_hg38_v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
      --knownSites /Ref_hg38_v0/Homo_sapiens_assembly38.known_indels.vcf.gz \
      --knownSites /Ref_hg38_v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
      --interval /Ref_hg38_v0/${INTERVAL} \
      --read-group-sm ${sample_id} \
      --read-group-lb ${sample_id} \
      --read-group-pl Illumina \
      --read-group-id-prefix ${sample_id} \
      --out-bam /Workdir/${sample_id}_dedup_sorted.bam \
      --out-qc-metrics-dir /Workdir/${sample_id}-qc-metrics/ \
      --out-recal-file /Workdir/${sample_id}-recal.txt \
      --tmp-dir /Workdir

  # ApplyBQSR
  apptainer exec --nv \
    -B /common/db/human_ref/hg38/parabricks:/Ref_hg38 \
    -B ${output_dir}:/InputDir \
    -B ${output_dir}:/OutputDir \
    -B /common/db/human_ref/hg38/v0:/Ref_hg38_v0 \
    /common/sif/clara-parabricks/4.4.0.sif pbrun \
      applybqsr \
      --ref /Ref_hg38/Homo_sapiens_assembly38.fasta \
      --in-bam /InputDir/${sample_id}_dedup_sorted.bam \
      --in-recal-file /InputDir/${sample_id}-recal.txt  \
      --interval /Ref_hg38_v0/${INTERVAL} \
      --out-bam /OutputDir/${sample_id}_recal.bam \
      --tmp-dir /OutputDir

  # Calculate Flagstats
  samtools flagstats -@ "${THREADS}" ${output_dir}/${sample_id}_recal.bam > ${output_dir}/${sample_id}_flagstat.txt

  # Convert BAM to CRAM
  samtools view -@ "${THREADS}" -C -T ${FASTA_FILE} \
    -o ${output_dir}/${sample_id}_dedup_sorted.cram ${output_dir}/${sample_id}_dedup_sorted.bam
  #samtools view -@ "${THREADS}" -C -T ${FASTA_FILE} \
  #  -o ${output_dir}/${sample_id}_recal.cram ${output_dir}/${sample_id}_recal.bam

  samtools index -@ "${THREADS}" -b ${output_dir}/${sample_id}_dedup_sorted.cram
  #samtools index -@ "${THREADS}" -b ${output_dir}/${sample_id}_recal.cram

  rm ${output_dir}/${sample_id}_dedup_sorted.bam* #${output_dir}/${sample_id}_recal.bam*
}

# --- Function: Run MosDepth ---
run_mosdepth() {
  ml purge
  ml apptainer

  local sample_id=$1
  local input_dir=$2
  local output_dir="${3}/${sample_id}"

  echo "LOG: Running MosDepth for sample ${sample_id}..."
  mkdir -p "${output_dir}" "${REPORT_DIR}"

  if [[ $SEQ_PLATFORM == "WGS" ]]
  then
    # MosDepth For WGS
    apptainer run --nv \
      -B ${output_dir}:/InputDir \
      -B ${output_dir}:/OutputDir \
      /common/sif/mosdepth/Mosdepth_0.3.10--h4e814b3_1.sif \
      mosdepth -n --fast-mode -t "${THREADS}" \
      --by 1000 \
      /OutputDir/${sample_id}_WGS_depth /InputDir/${sample_id}_recal.bam
  elif [[ $SEQ_PLATFORM == "WES" ]]
  then
    # MosDepth For WES or Target
    apptainer run --nv \
      -B ${output_dir}:/InputDir \
      -B ${output_dir}:/OutputDir \
      -B /project/o260003_CRU_BI/human_ref:/Ref \
      /common/sif/mosdepth/Mosdepth_0.3.10--h4e814b3_1.sif \
      mosdepth -n --fast-mode -t "${THREADS}" \
      --by /Ref/exome_calling_regions.v1.interval_list.bed.gz \
      /OutputDir/${sample_id}_WGS_depth /InputDir/${sample_id}_recal.bam
  fi
}

run_multiqc() {
  ml purge
  ml multiqc

  echo "LOG: Running MultiQC..."
  multiqc ${OUT_DIR} -o ${OUT_DIR}

  python ${QC_SCRIPT} ${OUT_DIR}/multiqc_data/multiqc_data.json ${OUT_DIR}/QC_SUMMARY.xlsx
}


# 3. RUNNING (Main Logic)
main() {
  # New directories
  CLEAN_DIR="${OUT_DIR}/01_CLEAN_QC"
  PREPROCESS_DIR="${OUT_DIR}/02_PREPROCESS"

  # Ensure directories exist
  mkdir -p "${OUT_DIR}" "${REPORT_DIR}" "${CLEAN_DIR}" "${PREPROCESS_DIR}"

  # Loop through read sample_sheet.csv file
  {
    # Skip header
    read -r header 

    # Extract 3 columns: ID, R1, and R2
    while IFS=',' read -r col1 col2 col3 || [ -n "$col1" ]; do
      
      base_name="${col1}"
      IN_DIR=$(dirname "$col2")

      if [[ "$base_name" =~ "WGS" ]]
      then
        SEQ_PLATFORM="WGS"
        INTERVAL="wgs_calling_regions.hg38.interval_list"
      elif [[ "$base_name" =~ "WES" ]]
      then
        SEQ_PLATFORM="WES"
        INTERVAL="exome_calling_regions.v1.interval_list"
      else
        echo "Wrong platform type for sample ${base_name}!!"
        exit
      fi
      
      local r1_file="${IN_DIR}/${base_name}_R1.fastq.gz"
      local r2_file="${IN_DIR}/${base_name}_R2.fastq.gz"

      # 0. Rename R1 & R2
      if [[ ! -f "$r1_file" ]]
      then
        mv "$col2" "$r1_file"
      fi
      
      if [[ ! -f "$r2_file" ]]
      then
        mv "$col3" "$r2_file"
      fi

      # 1. Initial QC
      run_fastqc "${r1_file}" "${r2_file}" "${REPORT_DIR}/fastqc_raw"

      # 2. Processing with fastp
      run_fastp "${r1_file}" "${r2_file}" "${base_name}" "${CLEAN_DIR}"

      # 3. Post-processing QC
      run_fastqc "${CLEAN_DIR}/${base_name}_R1_clean.fastq.gz" \
        "${CLEAN_DIR}/${base_name}_R2_clean.fastq.gz" \
        "${REPORT_DIR}/fastqc_clean"

      # 4. Preprocessing with FQ2BAM
      run_fq2bam "${base_name}" "${CLEAN_DIR}" "${PREPROCESS_DIR}"

      # 5. Calculate Coverage
      run_mosdepth "${base_name}" "${CLEAN_DIR}" "${PREPROCESS_DIR}"

    done

  } < "$SAMPLES_SHEET"

  run_multiqc

  echo "SUCCESS: Pipeline completed."
}

# Execution
main "$@"
