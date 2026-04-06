#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=47
#SBATCH --mem=150G
#SBATCH --gpus=1

# ==============================================================================
# Script Name : 02_run_mutect2.sh
# Description : Run Parabricks Mutect2 on one tumor-normal pair
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
readonly VCF_DIR="${CASE_DIR}/vcf"
readonly TMP_DIR="${CASE_DIR}/tmp_mutect2"

readonly REF_DIR="/common/db/human_ref/hg38/parabricks"
readonly PARABRICKS_SIF="/common/sif/clara-parabricks/4.5.1-1.sif"
readonly REF="${REF_DIR}/Homo_sapiens_assembly38.fasta"

# Assign BED here
readonly BED="/project/o260003_CRU_BI/human_ref/S33266436_Padded.chr.clean.bed"

readonly UNFILTERED_VCF="${VCF_DIR}/${CASE_ID}.unfiltered.vcf.gz"


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
    [[ -f "${BED}" ]] || { echo "ERROR: BED file not found: ${BED}"; exit 1; }

    mkdir -p "${VCF_DIR}" "${TMP_DIR}" logs/mutect2
}

load_modules() {
    ml purge
    ml apptainer
    ml bcftools
}

run_mutect2() {
    echo "=================================================="
    echo "STEP 2: MUTECT2 CALLING"
    echo "Case ID          : ${CASE_ID}"
    echo "Tumor sample ID  : ${TUMOR_SAMPLE_ID}"
    echo "Normal sample ID : ${NORMAL_SAMPLE_ID}"
    echo "=================================================="

    echo ">>> Running Mutect2 for ${CASE_ID}"
    apptainer exec --nv \
        -B "${REF_DIR},${OUT_DIR},$(dirname "${TUMOR_BAM}"),$(dirname "${NORMAL_BAM}"),$(dirname "${BED}")" \
        "${PARABRICKS_SIF}" \
        pbrun mutectcaller \
        --ref "${REF}" \
        --in-tumor-bam "${TUMOR_BAM}" \
        --in-normal-bam "${NORMAL_BAM}" \
        --interval "${BED}" \
        --tumor-name "${TUMOR_SAMPLE_ID}" \
        --normal-name "${NORMAL_SAMPLE_ID}" \
        --out-vcf "${UNFILTERED_VCF}" \
        --tmp-dir "${TMP_DIR}"

    echo ">>> Variant calling complete."

    if [[ -f "${UNFILTERED_VCF}" ]]; then
        local unfiltered_variants
        unfiltered_variants=$(bcftools view -H "${UNFILTERED_VCF}" | wc -l)
        echo ">>> ${unfiltered_variants} unfiltered variants found in ${UNFILTERED_VCF}"
    else
        echo "ERROR: Unfiltered VCF was not generated."
        exit 1
    fi
}

main() {
    check_inputs
    load_modules
    run_mutect2
}

main "$@"