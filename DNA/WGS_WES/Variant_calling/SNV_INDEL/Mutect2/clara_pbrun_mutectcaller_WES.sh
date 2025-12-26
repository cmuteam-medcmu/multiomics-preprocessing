#!/bin/bash

#SBATCH --account=o250027              
#SBATCH --job-name=mutectcaller_WES
#SBATCH --partition=gpu
#SBATCH --time=24:00:00 
#SBATCH --cpus-per-task=47
#SBATCH --mem=240G
#SBATCH --output=logs/mutect2/mutectcaller_WES_%j.log          
#SBATCH --error=logs/mutect2/mutectcaller_WES_%j.err
#SBATCH --gpus=1

ml apptainer
ml bcftools

# ===== User-defined Variables =====
PATIENT=$1
# The PATIENT variable must be set before running this script.
# Example: export PATIENT="patient_sample_name"
if [ -z "$PATIENT" ]; then
    echo "Error: The PATIENT environment variable is not set."
    exit 1
fi

# ===== Centralized Directories =====
PROJECT_DIR=/path_to_project_DIR
REF_DIR=/path_to_reference_DIR                  #/common/db/human_ref/hg38/parabricks
BAM_DIR=${PROJECT_DIR}/path_to_BAM_storage
OUT_DIR=${PROJECT_DIR}/Somatic_callers/${PATIENT}/vcf/mutect2

mkdir -p ${OUT_DIR} logs/mutect2

# ===== Input BAM Files =====
BAM_TUMOR=${BAM_DIR}/${PATIENT}_tumor_WES.recal.bam
BAM_NORMAL=${BAM_DIR}/${PATIENT}_normal_WES.recal.bam
# BAM_NORMAL=${BAM_DIR}/${PATIENT}_normal_WGS.exome.recal.bam

# ===== Reference and BED =====
REF=${REF_DIR}/Homo_sapiens_assembly38.fasta
BED=/path_to_Exon_BED_file.bed                #paht your exon captured bed file

# ===== Tumor / Normal Name =====
tumor=${PATIENT}_Tumor
normal=${PATIENT}_Normal

# ===== Output =====
UNFILTERED_VCF=${OUT_DIR}/${PATIENT}.unfiltered.vcf.gz

# ===== Run Mutect2 =====
echo ">>> Running Mutect2 for ${PATIENT}"

apptainer exec --nv \
  -B ${REF_DIR} \
  -B ${PROJECT_DIR} \
  /common/sif/clara-parabricks/4.5.1-1.sif pbrun \
  mutectcaller \
  --ref ${REF} \
  --in-tumor-bam ${BAM_TUMOR} \
  --in-normal-bam ${BAM_NORMAL} \
  --interval ${BED} \
  --tumor-name ${tumor} \
  --normal-name ${normal} \
  --out-vcf ${UNFILTERED_VCF}

echo ">>> Variant calling complete."

# Check if the VCF was created before trying to count variants
if [ -f "${UNFILTERED_VCF}" ]; then
    UNFILTERED_VARIANTS=$(bcftools view -H ${UNFILTERED_VCF} | wc -l)
    echo ">>> ${UNFILTERED_VARIANTS} unfiltered variants found in ${UNFILTERED_VCF}"
else
    echo "❌ ERROR: Unfiltered VCF was not generated. Check the error logs."
    exit 1
fi
