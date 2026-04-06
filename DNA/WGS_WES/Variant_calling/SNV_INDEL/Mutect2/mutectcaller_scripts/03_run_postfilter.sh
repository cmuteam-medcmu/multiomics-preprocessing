#!/bin/bash

#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# ==============================================================================
# Script Name : 03_run_postfilter.sh
# Description : Run FilterMutectCalls and optional custom somatic QC filtering
# ==============================================================================

set -euo pipefail

# =========================
# 1. INPUT ARGUMENTS
# =========================
readonly CASE_ID="${1:-}"
readonly TUMOR_SAMPLE_ID="${2:-}"
readonly NORMAL_SAMPLE_ID="${3:-}"
readonly FILTER_MODE="${4:-gatk_pass}"
# allowed: gatk_pass, custom_qc, all

# =========================
# 2. VARIABLE SETTING
# =========================
readonly OUT_DIR="/project/o260003_CRU_BI/MUTECT2_test"
readonly CASE_DIR="${OUT_DIR}/${CASE_ID}"
readonly RECAL_DIR="${CASE_DIR}/recal_tables"
readonly VCF_DIR="${CASE_DIR}/vcf"

readonly REF_DIR="/common/db/human_ref/hg38/parabricks"
readonly GATK_SIF="/common/sif/gatk/4.6.1.0.sif"
readonly REF="${REF_DIR}/Homo_sapiens_assembly38.fasta"

readonly UNFILTERED_VCF="${VCF_DIR}/${CASE_ID}.unfiltered.vcf.gz"
readonly CONTAMINATION_TABLE="${RECAL_DIR}/${CASE_ID}.contamination.table"
readonly SEGMENTS_TABLE="${RECAL_DIR}/${CASE_ID}.segments.table"

# Common GATK outputs
readonly FILTERED_VCF="${VCF_DIR}/${CASE_ID}.filtered.vcf.gz"
readonly GATK_PASS_VCF="${VCF_DIR}/${CASE_ID}.filtered.PASS.vcf.gz"

# Optional custom QC outputs
readonly FLAGGED_VCF="${VCF_DIR}/${CASE_ID}.quality_flagged.vcf.gz"
readonly FINAL_QUALITY_VCF="${VCF_DIR}/${CASE_ID}.filtered.PASS.QUALITY.vcf.gz"

# Custom thresholds
readonly TUMOR_DP_MIN=20
readonly NORMAL_DP_MIN=10
readonly TUMOR_ALT_AD_MIN=5

# =========================
# 3. FUNCTIONS
# =========================
check_inputs() {
    [[ -n "${CASE_ID}" ]] || { echo "ERROR: CASE_ID is required."; exit 1; }
    [[ -n "${TUMOR_SAMPLE_ID}" ]] || { echo "ERROR: TUMOR_SAMPLE_ID is required."; exit 1; }
    [[ -n "${NORMAL_SAMPLE_ID}" ]] || { echo "ERROR: NORMAL_SAMPLE_ID is required."; exit 1; }

    case "${FILTER_MODE}" in
        gatk_pass|custom_qc|all) ;;
        *)
            echo "ERROR: Invalid FILTER_MODE '${FILTER_MODE}'"
            echo "Allowed values: gatk_pass, custom_qc, all"
            exit 1
            ;;
    esac

    [[ -f "${REF}" ]] || { echo "ERROR: Reference FASTA not found: ${REF}"; exit 1; }
    [[ -f "${UNFILTERED_VCF}" ]] || { echo "ERROR: Unfiltered VCF not found: ${UNFILTERED_VCF}"; exit 1; }
    [[ -f "${CONTAMINATION_TABLE}" ]] || { echo "ERROR: Contamination table not found: ${CONTAMINATION_TABLE}"; exit 1; }
    [[ -f "${SEGMENTS_TABLE}" ]] || { echo "ERROR: Segments table not found: ${SEGMENTS_TABLE}"; exit 1; }

    mkdir -p logs/mutect2/gatk_filter logs/mutect2/somatic_qc
}

load_modules() {
    ml purge
    ml apptainer
    ml bcftools
    ml htslib
}

run_gatk_filter() {
    echo "=================================================="
    echo "STEP 3A: GATK FILTER"
    echo "Case ID     : ${CASE_ID}"
    echo "Filter mode : ${FILTER_MODE}"
    echo "=================================================="

    echo ">>> Running FilterMutectCalls..."
    apptainer exec \
        -B "${REF_DIR},${OUT_DIR}" \
        "${GATK_SIF}" \
        gatk FilterMutectCalls \
        -R "${REF}" \
        -V "${UNFILTERED_VCF}" \
        --contamination-table "${CONTAMINATION_TABLE}" \
        --tumor-segmentation "${SEGMENTS_TABLE}" \
        -O "${FILTERED_VCF}"

    echo ">>> GATK filtered VCF: ${FILTERED_VCF}"

    local filtered_count
    filtered_count=$(bcftools view -H "${FILTERED_VCF}" | wc -l)
    echo ">>> Variants in filtered VCF: ${filtered_count}"
}

extract_gatk_pass() {
    echo "=================================================="
    echo "STEP 3B: EXTRACT GATK PASS"
    echo "Case ID : ${CASE_ID}"
    echo "=================================================="

    [[ -f "${FILTERED_VCF}" ]] || { echo "ERROR: Filtered VCF not found: ${FILTERED_VCF}"; exit 1; }

    bcftools view -f PASS "${FILTERED_VCF}" -Oz -o "${GATK_PASS_VCF}"
    tabix -f -p vcf "${GATK_PASS_VCF}"

    local pass_count
    pass_count=$(bcftools view -H "${GATK_PASS_VCF}" | wc -l)

    echo ">>> GATK PASS VCF: ${GATK_PASS_VCF}"
    echo ">>> Final PASS variant count: ${pass_count}"
}

detect_sample_indices() {
    mapfile -t samples < <(bcftools query -l "${FILTERED_VCF}")

    if [[ "${#samples[@]}" -ne 2 ]]; then
        echo "ERROR: Expected exactly 2 samples in VCF, found ${#samples[@]}"
        printf '%s\n' "${samples[@]}"
        exit 1
    fi

    echo ">>> VCF sample order:"
    for i in "${!samples[@]}"; do
        echo ">>> [${i}] ${samples[$i]}"
    done

    TUMOR_IDX=""
    NORMAL_IDX=""

    for i in "${!samples[@]}"; do
        if [[ "${samples[$i]}" == "${TUMOR_SAMPLE_ID}" ]]; then
            TUMOR_IDX="${i}"
        elif [[ "${samples[$i]}" == "${NORMAL_SAMPLE_ID}" ]]; then
            NORMAL_IDX="${i}"
        fi
    done

    if [[ -z "${TUMOR_IDX}" || -z "${NORMAL_IDX}" ]]; then
        echo "ERROR: Could not match tumor/normal sample IDs in VCF."
        echo "Expected tumor  : ${TUMOR_SAMPLE_ID}"
        echo "Expected normal : ${NORMAL_SAMPLE_ID}"
        echo "VCF contains:"
        printf '%s\n' "${samples[@]}"
        exit 1
    fi

    echo ">>> Detected tumor sample index : [${TUMOR_IDX}] ${TUMOR_SAMPLE_ID}"
    echo ">>> Detected normal sample index: [${NORMAL_IDX}] ${NORMAL_SAMPLE_ID}"
}

run_custom_qc() {
    echo "=================================================="
    echo "STEP 3C: OPTIONAL CUSTOM QC FILTER"
    echo "Case ID          : ${CASE_ID}"
    echo "Tumor sample ID  : ${TUMOR_SAMPLE_ID}"
    echo "Normal sample ID : ${NORMAL_SAMPLE_ID}"
    echo "=================================================="

    [[ -f "${FILTERED_VCF}" ]] || { echo "ERROR: Filtered VCF not found: ${FILTERED_VCF}"; exit 1; }

    detect_sample_indices

    echo ">>> Applying optional custom filters..."
    echo ">>> Thresholds:"
    echo "    Tumor DP      >= ${TUMOR_DP_MIN}"
    echo "    Normal DP     >= ${NORMAL_DP_MIN}"
    echo "    Tumor alt AD  >= ${TUMOR_ALT_AD_MIN}"

    bcftools filter -m + -s LOW_T_DP -e "FORMAT/DP[${TUMOR_IDX}] < ${TUMOR_DP_MIN}" -Ou "${FILTERED_VCF}" \
    | bcftools filter -m + -s LOW_N_DP -e "FORMAT/DP[${NORMAL_IDX}] < ${NORMAL_DP_MIN}" -Ou \
    | bcftools filter -m + -s LOW_T_AD -e "FORMAT/AD[${TUMOR_IDX}:1] < ${TUMOR_ALT_AD_MIN}" -Oz -o "${FLAGGED_VCF}"

    tabix -f -p vcf "${FLAGGED_VCF}"

    echo ""
    echo ">>> FILTER SUMMARY (Audit Trail):"
    echo "------------------------------------------"
    echo "Count | Filter Status"
    bcftools query -f '%FILTER\n' "${FLAGGED_VCF}" | sort | uniq -c
    echo "------------------------------------------"

    echo ">>> Extracting final PASS variants after optional custom QC..."
    bcftools view -f PASS "${FLAGGED_VCF}" -Oz -o "${FINAL_QUALITY_VCF}"
    tabix -f -p vcf "${FINAL_QUALITY_VCF}"

    local final_count
    final_count=$(bcftools view -H "${FINAL_QUALITY_VCF}" | wc -l)

    echo ""
    echo ">>> Custom QC complete."
    echo ">>> Flagged VCF: ${FLAGGED_VCF}"
    echo ">>> Final High-Confidence VCF: ${FINAL_QUALITY_VCF}"
    echo ">>> Final Somatic Variant Count: ${final_count}"
}

main() {
    check_inputs
    load_modules
    run_gatk_filter

    case "${FILTER_MODE}" in
        gatk_pass)
            extract_gatk_pass
            ;;
        custom_qc)
            run_custom_qc
            ;;
        all)
            extract_gatk_pass
            run_custom_qc
            ;;
    esac
}

main "$@"