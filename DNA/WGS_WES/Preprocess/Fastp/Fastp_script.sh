#!/bin/bash

#SBATCH --account=o250022               
#SBATCH --job-name=Fastp 
#SBATCH --partition=compute
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=10G
#SBATCH --array=0-11%6
#SBATCH --output=Fastp_output_%A_%a.log
#SBATCH --error=Fastp_error_%A_%a.log

ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" /project/o250003_CFBile/Target_Seq/02_Macrogen_Prep/00_RAW_FASTQ/Round1_id.txt)

INPUT_DIR='/project/o250003_CFBile/Target_Seq/02_Macrogen_Prep/00_RAW_FASTQ'
OUTPUT_DIR='/project/o250003_CFBile/Target_Seq/02_Macrogen_Prep/01_CLEAN_QC'

ml purge
ml cutadapt
ml fastp

## For Defalt
fastp -i ${INPUT_DIR}/${ID}_1.fastq.gz \
    -I ${INPUT_DIR}/${ID}_2.fastq.gz \
    -o ${OUTPUT_DIR}/${ID}_1_cleaned.fastq.gz \
    -O ${OUTPUT_DIR}/${ID}_2_cleaned.fastq.gz \
    -w 20 --trim_poly_g \
    --html ${OUTPUT_DIR}/${ID}_fastp.html \
    --json ${OUTPUT_DIR}/${ID}_fastp.json \
    --dont_overwrite

## For UMI
# fastp -i ${INPUT_DIR}/${ID}_1.fastq.gz \
#     -I ${INPUT_DIR}/${ID}_2.fastq.gz \
#     -o ${OUTPUT_DIR}/${ID}_1_cleaned.fastq.gz \
#     -O ${OUTPUT_DIR}/${ID}_2_cleaned.fastq.gz \
#     -U --umi_loc=per_read --umi_len=8 \
#     -w 20 --trim_poly_g \
#     --html ${OUTPUT_DIR}/${ID}_fastp.html \
#     --json ${OUTPUT_DIR}/${ID}_fastp.json \
#     --dont_overwrite

ml purge
ml fastqc

fastqc ${OUTPUT_DIR}/${ID}_1_cleaned.fastq.gz ${OUTPUT_DIR}/${ID}_2_cleaned.fastq.gz -t 2 --outdir ${OUTPUT_DIR}

# ml purge
# ml multiqc

# multiqc ${OUTPUT_DIR}/
