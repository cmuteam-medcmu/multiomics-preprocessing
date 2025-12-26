#!/bin/bash

#SBATCH --account=o250027              
#SBATCH --job-name=pbrun_somatic
#SBATCH --partition=gpu
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --output=logs/mutect2/pbrun_somatic_%j.log          
#SBATCH --error=logs/mutect2/pbrun_somatic_%j.err
#SBATCH --gpus=1

# Load required modules
ml apptainer
ml bcftools

# ===== Patient ID =====
PATIENT=$1

if [ -z "${PATIENT}" ]; then
    echo "ERROR: No patient ID provided. Usage: $0 <LOT> <PATIENT_ID>"
    exit 1
fi

# ===== Centralized Directories =====
PROJECT_DIR=/path_to_project_DIR
REF_DIR=/path_to_reference_DIR                  #/common/db/human_ref/hg38/parabricks
KNOWN_SITES_DIR=/path_to_KNOWN_SITES_DIR        #/common/db/human_ref/hg38/v0
FASTQ_DIR=${PROJECT_DIR}/path_to_fp.fastq_DIR
BAM_DIR=${PROJECT_DIR}/path_to_BAM_storage_DIR

# A single main output directory for this patient's somatic analysis
OUT_DIR=${PROJECT_DIR}/Somatic_callers/${PATIENT}
OUT_VCF_DIR=${OUT_DIR}/vcf/mutect2
OUT_RECAL_DIR=${OUT_DIR}/recal_tables
TMP_DIR=${OUT_DIR}/tmp

# Create all necessary output directories at once
mkdir -p ${BAM_DIR} ${OUT_DIR} ${OUT_VCF_DIR} ${OUT_RECAL_DIR} ${TMP_DIR} logs

# ===== Reference and Annotation Files =====
REF=${REF_DIR}/Homo_sapiens_assembly38.fasta
KNOWN_SITES_SNP=${KNOWN_SITES_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz
KNOWN_SITES_INDEL=${KNOWN_SITES_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
BED=/path_to_Exon_BED_file.bed                #paht your exon captured bed file

# ===== Input FastQ Files =====
#Tumor fastq (WES)
FASTQ_T1=${FASTQ_DIR}/${PATIENT}_read_1.fp.fastq.gz
FASTQ_T2=${FASTQ_DIR}/${PATIENT}_read_2.fp.fastq.gz
#Normal fastq (WGS)
FASTQ_N1=${FASTQ_DIR}/${PATIENT}_read_1.fp.fastq.gz
FASTQ_N2=${FASTQ_DIR}/${PATIENT}_read_2.fp.fastq.gz

# ===== Dynamic Sample Names and Read Groups =====
TUMOR_NAME="${PATIENT}_Tumor"
NORMAL_NAME="${PATIENT}_Normal"
RG_TUMOR="@RG\tID:${TUMOR_NAME}\tSM:${TUMOR_NAME}\tPL:ILLUMINA\tLB:lib1\tPU:${TUMOR_NAME}.L1"
RG_NORMAL="@RG\tID:${NORMAL_NAME}\tSM:${NORMAL_NAME}\tPL:ILLUMINA\tLB:lib1\tPU:${NORMAL_NAME}.L1"

# ===== Output Files =====
BAM_TUMOR=${BAM_DIR}/${PATIENT}_tumor_WES.recal.bam
BAM_NORMAL=${BAM_DIR}/${PATIENT}_normal_WGS.recal.exome.bam
UNFILTERED_VCF=${OUT_VCF_DIR}/${PATIENT}.unfiltered.vcf.gz
TUMOR_RECAL_TXT=${OUT_RECAL_DIR}/${PATIENT}_tumor_WES.recal.txt
NORMAL_RECAL_TXT=${OUT_RECAL_DIR}/${PATIENT}_normal_WGS.recal.txt


# ===== Run Parabricks Somatic Pipeline =====
echo ">>> Starting Parabricks Somatic pipeline for patient ${PATIENT} (WGS vs WES)"

apptainer exec --nv \
  -B ${PROJECT_DIR} \
  -B ${REF_DIR} \
  /common/sif/clara-parabricks/4.5.1-1.sif pbrun \
  somatic \
    --ref ${REF} \
    --knownSites ${KNOWN_SITES_SNP} \
    --knownSites ${KNOWN_SITES_INDEL} \
    --interval ${BED} \
    --in-tumor-fq ${FASTQ_T1} ${FASTQ_T2} "${RG_TUMOR}" \
    --in-normal-fq ${FASTQ_N1} ${FASTQ_N2} "${RG_NORMAL}" \
    --bwa-options="-Y" \
    --out-vcf ${UNFILTERED_VCF} \
    --out-tumor-bam ${BAM_TUMOR} \
    --out-tumor-recal-file ${TUMOR_RECAL_TXT} \
    --out-normal-bam ${BAM_NORMAL} \
    --out-normal-recal-file ${NORMAL_RECAL_TXT} \
    --tmp-dir ${TMP_DIR}

echo ">>> Parabricks Somatic pipeline finished successfully."
rm -r ${TMP_DIR}

# Check if the VCF was created before trying to count variants
if [ -f "${UNFILTERED_VCF}" ]; then
    UNFILTERED_VARIANTS=$(bcftools view -H ${UNFILTERED_VCF} | wc -l)
    echo ">>> ${UNFILTERED_VARIANTS} unfiltered variants found in ${UNFILTERED_VCF}"
else
    echo "❌ ERROR: Unfiltered VCF was not generated. Check the error logs."
    exit 1
fi
