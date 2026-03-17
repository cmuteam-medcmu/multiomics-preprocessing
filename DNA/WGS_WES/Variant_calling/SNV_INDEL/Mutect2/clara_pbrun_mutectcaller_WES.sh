#!/bin/bash

#SBATCH --account=o250027              #Add your projectID
#SBATCH --job-name=pbrun_mutectcaller
#SBATCH --partition=gpu
#SBATCH --time=24:00:00 
#SBATCH --cpus-per-task=47
#SBATCH --mem=240G
#SBATCH --output=logs/mutect2/mutectcaller_WES_%j.log          
#SBATCH --error=logs/mutect2/mutectcaller_WES_%j.err
#SBATCH --gpus=1

ml apptainer
ml bcftools

# ===== Patient ID =====
PATIENT=$1

if [ -z "${PATIENT}" ]; then
    echo "ERROR: No patient ID provided."
    echo "Usage: sbatch $0 <PATIENT_ID>"
    exit 1
fi

# ===== Centralized Directories =====
PROJECT_DIR=/path_to_project_DIR
REF_DIR=/path_to_reference_DIR                      # e.g. /common/db/human_ref/hg38/parabricks
BAM_DIR=${PROJECT_DIR}/path_to_BAM_storage

OUT_DIR=${PROJECT_DIR}/Somatic_callers/${PATIENT}
OUT_VCF_DIR=${OUT_DIR}/vcf/mutect2
TMP_DIR=${OUT_DIR}/tmp

mkdir -p ${OUT_VCF_DIR} ${TMP_DIR} logs/mutect2

# ===== Input BAM Files =====
# Use the WES-WES normal BAM by default.
# For a WGS Normal vs WES Tumor run, uncomment the WGS line and comment out the WES line.
BAM_TUMOR=${BAM_DIR}/${PATIENT}_tumor_WES.recal.bam
BAM_NORMAL=${BAM_DIR}/${PATIENT}_normal_WES.recal.bam
# BAM_NORMAL=${BAM_DIR}/${PATIENT}_normal_WGS.recal.bam   # WGS Normal (use with WGS-WES somatic output)

# ===== Reference and BED =====
REF=${REF_DIR}/Homo_sapiens_assembly38.fasta
BED=/path_to_Exon_BED_file.bed                      # path to your exon capture BED file

# ===== Tumor / Normal Sample Names =====
# These must match the SM: tag in the read group of the input BAMs
TUMOR_NAME="${PATIENT}_Tumor"
NORMAL_NAME="${PATIENT}_Normal"

# ===== Output =====
UNFILTERED_VCF=${OUT_VCF_DIR}/${PATIENT}.unfiltered.vcf.gz

# ===== Validate Inputs =====
for f in "${BAM_TUMOR}" "${BAM_NORMAL}" "${REF}" "${BED}"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required input file not found: $f"
        exit 1
    fi
done

# ===== Run Mutect2 =====
echo ">>> Running Mutect2 (mutectcaller) for ${PATIENT}"

apptainer exec --nv \
  -B ${REF_DIR} \
  -B ${PROJECT_DIR} \
  /common/sif/clara-parabricks/4.5.1-1.sif pbrun \
  mutectcaller \
    --ref ${REF} \
    --in-tumor-bam ${BAM_TUMOR} \
    --in-normal-bam ${BAM_NORMAL} \
    --interval ${BED} \
    --tumor-name ${TUMOR_NAME} \
    --normal-name ${NORMAL_NAME} \
    --out-vcf ${UNFILTERED_VCF} \
    --tmp-dir ${TMP_DIR}

PBRUN_EXIT=$?
if [ ${PBRUN_EXIT} -ne 0 ]; then
    echo "ERROR: mutectcaller failed with exit code ${PBRUN_EXIT}."
    exit ${PBRUN_EXIT}
fi

echo ">>> Variant calling complete."
rm -rf ${TMP_DIR}

# ===== Validate Output =====
if [ -f "${UNFILTERED_VCF}" ]; then
    UNFILTERED_VARIANTS=$(bcftools view -H ${UNFILTERED_VCF} | wc -l)
    echo ">>> ${UNFILTERED_VARIANTS} unfiltered variants found in ${UNFILTERED_VCF}"
else
    echo "ERROR: Unfiltered VCF was not generated. Check the error logs."
    exit 1
fi
