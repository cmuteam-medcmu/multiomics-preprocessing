#!/bin/bash
#SBATCH --job-name=QC_1case
#SBATCH --partition=gpu  #(short/compute/gpu)
#SBATCH --time=24:00:00  #hh:mm:ss
#SBATCH --mem=240G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/qc/QC_%j.out
#SBATCH --error=logs/qc/QC_%j.err

module load fastqc
module load fastp

#--- Define Paths---#
INPUT_DIR=/path_to_fastq_DIR
OUT_DIR_fastqc1=/path_to_1st_fastQC_output_DIR
OUT_DIR_fastp=/path_to_fastp_output_DIR
OUT_DIR_fastqc2=/path_to_2nd_fastQC_output_DIR

#Input files
INPUT_FILE_1=${INPUT_DIR}/read1.fastq.gz
INPUT_FILE_2=${INPUT_DIR}/read2.fastq.gz
#Output files (define file's name here)
OUTPUT_FILE_1=${OUT_DIR_fastp}/read1.fp.fastq.gz
OUTPUT_FILE_2=${OUT_DIR_fastp}/read2.fp.fastq.gz
HTML_FILE=${OUT_DIR_fastp}/file.fastp.html
JSON_FILE=${OUT_DIR_fastp}/file.fastp.json

#--- Making OUT_DIRS ----#
mkdir -p ${OUT_DIR_fastqc1} ${OUT_DIR_fastqc2} ${OUT_DIR_fastp} logs/qc

#---- 1st QC ---#
fastqc -o ${OUT_DIR_fastqc1} ${INPUT_FILE_1}
fastqc -o ${OUT_DIR_fastqc1} ${INPUT_FILE_2}

#---- cut adapter & trim poly G by fastp ---#
fastp -i ${INPUT_FILE_1} -I ${INPUT_FILE_2} \
          -o ${OUTPUT_FILE_1} \
          -O ${OUTPUT_FILE_2} \
          -h ${HTML_FILE}\
          -j ${JSON_FILE}\
          --trim_poly_g

#---- 2nd QC ---#
#qc after fastp
fastqc -o ${OUT_DIR_fastqc2} ${OUTPUT_FILE_1}
fastqc -o ${OUT_DIR_fastqc2} ${OUTPUT_FILE_2}

echo "QC finised, fastq files ready to use"
