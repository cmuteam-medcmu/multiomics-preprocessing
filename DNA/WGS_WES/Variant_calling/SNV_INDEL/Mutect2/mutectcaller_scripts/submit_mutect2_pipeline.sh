#!/bin/bash

# ==============================================================================
# Script Name : submit_mutect2_pipeline.sh
# Description : Submit shared somatic pipeline by reading pairing sample sheet
# ==============================================================================

set -euo pipefail

# =========================
# 1. VARIABLE SETTING
# =========================
readonly SAMPLES_SHEET="/project/o260003_CRU_BI/MUTECT2_test/sample_sheet_multisamples.csv"
readonly SCRIPT_DIR="/project/o260003_CRU_BI/script/mutectcaller_scripts"
readonly LOG_ROOT="/project/o260003_CRU_BI/MUTECT2_test/logs"

readonly STEP1="${SCRIPT_DIR}/01_run_gatk_tables.sh"
readonly STEP2="${SCRIPT_DIR}/02_run_mutect2.sh"
readonly STEP3="${SCRIPT_DIR}/03_run_postfilter.sh"

readonly SLURM_ACCOUNT="o260003"

# allowed: gatk_pass, custom_qc, all
readonly FILTER_MODE="all"

# allowed: yes, no
# yes = next case waits until previous case finishes STEP3
# no  = all cases are submitted independently
readonly SERIAL_BY_CASE="yes"

# =========================
# 2. FUNCTIONS
# =========================
trim() {
    local var="$1"
    var="${var#"${var%%[![:space:]]*}"}"
    var="${var%"${var##*[![:space:]]}"}"
    printf '%s' "${var}"
}

check_inputs() {
    [[ -f "${SAMPLES_SHEET}" ]] || { echo "ERROR: Sample sheet not found: ${SAMPLES_SHEET}"; exit 1; }

    for f in "${STEP1}" "${STEP2}" "${STEP3}"; do
        [[ -f "${f}" ]] || { echo "ERROR: Script not found: ${f}"; exit 1; }
    done

    case "${FILTER_MODE}" in
        gatk_pass|custom_qc|all) ;;
        *)
            echo "ERROR: Invalid FILTER_MODE '${FILTER_MODE}'"
            echo "Allowed values: gatk_pass, custom_qc, all"
            exit 1
            ;;
    esac

    case "${SERIAL_BY_CASE}" in
        yes|no) ;;
        *)
            echo "ERROR: Invalid SERIAL_BY_CASE '${SERIAL_BY_CASE}'"
            echo "Allowed values: yes, no"
            exit 1
            ;;
    esac

    mkdir -p "${LOG_ROOT}"
}

submit_pipeline() {
    local case_id=$1
    local tumor_sample_id=$2
    local normal_sample_id=$3
    local tumor_bam=$4
    local normal_bam=$5
    local previous_case_final_jid=${6:-}

    local jid1 jid2 jid3
    local log_dir_step1 log_dir_step2 log_dir_step3
    local dep_args=()

    if [[ "${SERIAL_BY_CASE}" == "yes" && -n "${previous_case_final_jid}" ]]; then
        dep_args=(--dependency="afterok:${previous_case_final_jid}")
    fi

    log_dir_step1="${LOG_ROOT}/${case_id}/gatk_table"
    mkdir -p "${log_dir_step1}"

    jid1=$(sbatch --parsable \
        --account="${SLURM_ACCOUNT}" \
        "${dep_args[@]}" \
        --job-name="gatk_${case_id}" \
        --output="${log_dir_step1}/gatk_tables_%j.log" \
        --error="${log_dir_step1}/gatk_tables_%j.err" \
        "${STEP1}" \
        "${case_id}" "${tumor_sample_id}" "${normal_sample_id}" "${tumor_bam}" "${normal_bam}")

    log_dir_step2="${LOG_ROOT}/${case_id}/mutect2"
    mkdir -p "${log_dir_step2}"

    jid2=$(sbatch --parsable \
        --account="${SLURM_ACCOUNT}" \
        "${dep_args[@]}" \
        --job-name="mutect_${case_id}" \
        --output="${log_dir_step2}/mutectcaller_%j.log" \
        --error="${log_dir_step2}/mutectcaller_%j.err" \
        "${STEP2}" \
        "${case_id}" "${tumor_sample_id}" "${normal_sample_id}" "${tumor_bam}" "${normal_bam}")

    log_dir_step3="${LOG_ROOT}/${case_id}/postfilter"
    mkdir -p "${log_dir_step3}"

    jid3=$(sbatch --parsable \
        --account="${SLURM_ACCOUNT}" \
        --dependency="afterok:${jid1}:${jid2}" \
        --job-name="post_${case_id}" \
        --output="${log_dir_step3}/postfilter_%j.log" \
        --error="${log_dir_step3}/postfilter_%j.err" \
        "${STEP3}" \
        "${case_id}" "${tumor_sample_id}" "${normal_sample_id}" "${FILTER_MODE}")

    echo "Case ID      : ${case_id}"
    echo "  GATK tables: ${jid1}"
    echo "  Mutect2    : ${jid2}"
    echo "  Postfilter : ${jid3}"
    echo "  Mode       : ${FILTER_MODE}"
    echo "  Serial     : ${SERIAL_BY_CASE}"
    if [[ "${SERIAL_BY_CASE}" == "yes" && -n "${previous_case_final_jid}" ]]; then
        echo "  Starts after previous case final job: ${previous_case_final_jid}"
    elif [[ "${SERIAL_BY_CASE}" == "yes" ]]; then
        echo "  Starts immediately (first case)"
    else
        echo "  Submitted independently"
    fi
    echo

    printf '%s\n' "${jid3}"
}

main() {
    check_inputs

    local previous_case_final_jid=""
    local current_case_final_jid=""

    {
        read -r header

        while IFS=',' read -r case_id tumor_sample_id normal_sample_id tumor_bam normal_bam || [[ -n "${case_id:-}" ]]; do
            case_id=$(trim "${case_id:-}")
            tumor_sample_id=$(trim "${tumor_sample_id:-}")
            normal_sample_id=$(trim "${normal_sample_id:-}")
            tumor_bam=$(trim "${tumor_bam:-}")
            normal_bam=$(trim "${normal_bam:-}")

            [[ -z "${case_id}" ]] && continue

            current_case_final_jid=$(submit_pipeline \
                "${case_id}" \
                "${tumor_sample_id}" \
                "${normal_sample_id}" \
                "${tumor_bam}" \
                "${normal_bam}" \
                "${previous_case_final_jid}" | tail -n 1)

            if [[ "${SERIAL_BY_CASE}" == "yes" ]]; then
                previous_case_final_jid="${current_case_final_jid}"
            fi
        done
    } < "${SAMPLES_SHEET}"
}

main "$@"