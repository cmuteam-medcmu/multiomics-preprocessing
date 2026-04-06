#!/bin/bash

#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# ==============================================================================
# Script Name : 01_run_gatk_tables.sh
# Description : Generate pileup and contamination tables for FilterMutectCalls
# ==============================================================================

set -euo pipefail

# =========================
# 1. INPUT ARGUMENTS
# =========================
readonly CASE_ID="${1:-}"
readonly TUMOR_SAMPLE_ID="${2:-}"
readonly NORMAL_SAMPLE_ID="${3:-}"
readonly TUMOR_BAM="${4:-}"
readonly NORMAL_BAM="${5:-}"

# =========================
# 2. VARIABLE SETTING
# =========================
readonly OUT_DIR="/project/o260003_CRU_BI/MUTECT2_test"
readonly CASE_DIR="${OUT_DIR}/${CASE_ID}"
readonly RECAL_DIR="${CASE_DIR}/recal_tables"

readonly REF_DIR="/common/db/human_ref/hg38/parabricks"
readonly KNOWN_SITES_DIR="/common/db/human_ref/hg38/v0"
readonly GATK_SIF="/common/sif/gatk/4.6.1.0.sif"

readonly REF="${REF_DIR}/Homo_sapiens_assembly38.fasta"
readonly KNOWN_SITES_SNP="${KNOWN_SITES_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"

# Assign BED here
readonly BED="/project/o260003_CRU_BI/human_ref/S33266436_Padded.chr.clean.bed"

# Outputs
readonly TUMOR_PILEUPS_TABLE="${RECAL_DIR}/${CASE_ID}.tumor.pileups.table"
readonly NORMAL_PILEUPS_TABLE="${RECAL_DIR}/${CASE_ID}.normal.pileups.table"
readonly CONTAMINATION_TABLE="${RECAL_DIR}/${CASE_ID}.contamination.table"
readonly SEGMENTS_TABLE="${RECAL_DIR}/${CASE_ID}.segments.table"

# =========================
# 3. FUNCTIONS
# =========================
check_inputs() {
    [[ -n "${CASE_ID}" ]] || { echo "ERROR: CASE_ID is required."; exit 1; }
    [[ -n "${TUMOR_SAMPLE_ID}" ]] || { echo "ERROR: TUMOR_SAMPLE_ID is required."; exit 1; }
    [[ -n "${NORMAL_SAMPLE_ID}" ]] || { echo "ERROR: NORMAL_SAMPLE_ID is required."; exit 1; }
    [[ -f "${TUMOR_BAM}" ]] || { echo "ERROR: Tumor BAM not found: ${TUMOR_BAM}"; exit 1; }
    [[ -f "${NORMAL_BAM}" ]] || { echo "ERROR: Normal BAM not found: ${NORMAL_BAM}"; exit 1; }
    [[ -f "${REF}" ]] || { echo "ERROR: Reference FASTA not found: ${REF}"; exit 1; }
    [[ -f "${KNOWN_SITES_SNP}" ]] || { echo "ERROR: Known sites VCF not found: ${KNOWN_SITES_SNP}"; exit 1; }
    [[ -f "${BED}" ]] || { echo "ERROR: BED file not found: ${BED}"; exit 1; }

    mkdir -p "${RECAL_DIR}" logs/mutect2/gatk_table
}

load_modules() {
    ml purge
    ml apptainer
}

run_gatk_tables() {
    echo "=================================================="
    echo "STEP 1: GATK TABLES"
    echo "Case ID          : ${CASE_ID}"
    echo "Tumor sample ID  : ${TUMOR_SAMPLE_ID}"
    echo "Normal sample ID : ${NORMAL_SAMPLE_ID}"
    echo "=================================================="

    echo ">>> Running GetPileupSummaries on tumor BAM..."
    apptainer exec \
        -B "${REF_DIR},${KNOWN_SITES_DIR},${OUT_DIR},$(dirname "${TUMOR_BAM}"),$(dirname "${NORMAL_BAM}"),$(dirname "${BED}")" \
        "${GATK_SIF}" \
        gatk GetPileupSummaries \
        -I "${TUMOR_BAM}" \
        -V "${KNOWN_SITES_SNP}" \
        -L "${BED}" \
        -R "${REF}" \
        -O "${TUMOR_PILEUPS_TABLE}"

    echo ">>> Running GetPileupSummaries on normal BAM..."
    apptainer exec \
        -B "${REF_DIR},${KNOWN_SITES_DIR},${OUT_DIR},$(dirname "${TUMOR_BAM}"),$(dirname "${NORMAL_BAM}"),$(dirname "${BED}")" \
        "${GATK_SIF}" \
        gatk GetPileupSummaries \
        -I "${NORMAL_BAM}" \
        -V "${KNOWN_SITES_SNP}" \
        -L "${BED}" \
        -R "${REF}" \
        -O "${NORMAL_PILEUPS_TABLE}"

    echo ">>> Running CalculateContamination..."
    apptainer exec \
        -B "${OUT_DIR}" \
        "${GATK_SIF}" \
        gatk CalculateContamination \
        -I "${TUMOR_PILEUPS_TABLE}" \
        -matched "${NORMAL_PILEUPS_TABLE}" \
        -O "${CONTAMINATION_TABLE}" \
        -segments "${SEGMENTS_TABLE}"

    echo ">>> GATK tables ready in ${RECAL_DIR}"
    echo ">>> Contamination table:"
    cat "${CONTAMINATION_TABLE}"
}

main() {
    check_inputs
    load_modules
    run_gatk_tables
}

main "$@"