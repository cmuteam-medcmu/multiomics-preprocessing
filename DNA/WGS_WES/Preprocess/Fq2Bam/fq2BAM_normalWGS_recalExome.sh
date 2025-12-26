#!/bin/bash

#SBATCH --account=o250027
#SBATCH --job-name=fq2bam_normal_WGS
#SBATCH --partition=gpu
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=240G
#SBATCH --output=logs/fq2bam_normal_WGS_%j.out
#SBATCH --error=logs/fq2bam_normal_WGS_%j.err
#SBATCH --gpus=1

PATIENT=$1

# Load required modules
ml purge
ml apptainer
ml samtools

# ===== User-defined Variables =====
# The PATIENT variable must be set before running this script.
# Example: export PATIENT="patient_sample_name"
if [ -z "$PATIENT" ]; then
    echo "Error: The PATIENT environment variable is not set."
    exit 1
fi

# ===== Input / Output Directories =====
PROJECT_DIR=/path_to_project_DIR
KNOWN_SITES_DIR=/path_to_KNOWN_SITES_DIR        #/common/db/human_ref/hg38/v0
REF_DIR=/path_to_reference_DIR                  #/common/db/human_ref/hg38/parabricks
FASTQ_DIR=${PROJECT_DIR}/path_to_fp.fastq_DIR
OUT_DIR=${PROJECT_DIR}/path_to_BAM_storage_DIR_or_output_DIR

# Create output directory if it doesn't exist
mkdir -p ${OUT_DIR} logs ${OUT_DIR}/tmp

# ===== Reference and Known Sites =====
FASTQ_1=${FASTQ_DIR}/read_1.fp.fastq.gz
FASTQ_2=${FASTQ_DIR}/read_2.fp.fastq.gz

KNOWN_SITES_SNP=${KNOWN_SITES_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz
KNOWN_SITES_INDEL=${KNOWN_SITES_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
REF=${REF_DIR}/Homo_sapiens_assembly38.fasta
EXON_BED=${PROJECT_DIR}/exon_capture_BED_file.bed

# ===== Output Files =====
OUTPUT_BAM=${OUT_DIR}/${PATIENT}_WGS.exome.bam
OUTPUT_RECAL_TXT=${OUT_DIR}/${PATIENT}_WGS.exome.txt
OUTPUT_EXOME_BAM=${OUT_DIR}/${PATIENT}_WGS.exome.recal.bam
OUTPUT_QC_DIR=${OUT_DIR}/${PATIENT}_WGS.exome-qc-metrics

# ===== Run fq2bam Pipeline =====
apptainer exec --nv \
  -B ${PROJECT_DIR} \
  -B ${KNOWN_SITES_DIR} \
  -B ${REF_DIR} \
  /common/sif/clara-parabricks/4.5.1-1.sif pbrun \
  fq2bam \
  --ref ${REF} \
  --in-fq ${FASTQ_1} ${FASTQ_2}  \
  --knownSites ${KNOWN_SITES_SNP} \
  --knownSites ${KNOWN_SITES_INDEL} \
  --interval-file ${EXON_BED} \
  --out-bam ${OUTPUT_BAM} \
  --out-recal-file ${OUTPUT_RECAL_TXT} \
  --out-qc-metrics-dir ${OUTPUT_QC_DIR} \
  --tmp-dir ${OUT_DIR}/tmp \
  --read-group-sm "${PATIENT}_Normal" \
  --read-group-id "${PATIENT}_Normal" \
  --read-group-lb "lib1" \
  --read-group-pl "ILLUMINA"

# ===== ApplyBQSR =====
apptainer exec --nv \
  -B ${REF_DIR} \
  -B ${PROJECT_DIR}\
  /common/sif/clara-parabricks/4.5.1-1.sif pbrun \
      applybqsr \
      --ref ${REF} \
      --in-bam ${OUTPUT_BAM} \
      --in-recal-file ${OUTPUT_RECAL_TXT}  \
      --out-bam ${OUTPUT_EXOME_BAM} \
      --tmp-dir ${OUT_DIR}/tmp

samtools index -@ 20 -b ${OUTPUT_EXOME_BAM}

rm -r ${OUT_DIR}/tmp

echo "Script finished successfully."
