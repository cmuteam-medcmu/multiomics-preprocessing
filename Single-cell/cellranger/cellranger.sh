#!/bin/bash

#SBATCH --job-name=Cellranger
#SBATCH --partition=compute
#SBATCH --time=05-00:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --output=/home/warunyoo.p/o250031_COVID68/cellranger_%j.out
#SBATCH --error=/home/warunyoo.p/o250031_COVID68/cellranger_%j.err
#SBATCH --account=o250031

#command: cellranger cloud auth setup
#register 10x account to get a cretdit token : https://cloud.10xgenomics.com/account/security
#Enter access token: 
#get access token path & copy to -tenx-cloud-token-path= in. command

module load cellranger/9.0.1

cellranger count --id=visit3-1 \
           --transcriptome=/common/db/human_ref/cellranger/refdata-gex-GRCh38-2024-A \
           --fastqs=/home/warunyoo.p/o250031_COVID68/HN00253241/HN00253241_10X_RawData_Outs/CMU_CV_B1-1_GEX/22YM3HLT3 \
           --sample=CMU_CV_B1-1_GEX \
           --create-bam=true \
           --localcores=64 \
           --localmem=256 \
           --tenx-cloud-token-path=/home/warunyoo.p/.config/txg/credentials

cellranger count --id=visit4-1 \
           --transcriptome=/common/db/human_ref/cellranger/refdata-gex-GRCh38-2024-A \
           --fastqs=/home/warunyoo.p/o250031_COVID68/HN00253241/HN00253241_10X_RawData_Outs/CMU_CV_B2-1_GEX/22YM3HLT3 \
           --sample=CMU_CV_B2-1_GEX \
           --create-bam=true \
           --localcores=64 \
           --localmem=256 \
           --tenx-cloud-token-path=/home/warunyoo.p/.config/txg/credentials

cellranger count --id=visit4-2 \
           --transcriptome=/common/db/human_ref/cellranger/refdata-gex-GRCh38-2024-A \
           --fastqs=/home/warunyoo.p/o250031_COVID68/HN00253241/HN00253241_10X_RawData_Outs/CMU_CV_B2-2_GEX/22YM3HLT3 \
           --sample=CMU_CV_B2-2_GEX \
           --create-bam=true \
           --localcores=64 \
           --localmem=256 \
           --tenx-cloud-token-path=/home/warunyoo.p/.config/txg/credentials

cellranger count --id=visit1-1 \
           --transcriptome=/common/db/human_ref/cellranger/refdata-gex-GRCh38-2024-A \
           --fastqs=/home/warunyoo.p/o250031_COVID68/HN00254955/HN00254955_10X_RawData_Outs/CMU_CV_B3-1_GEX/22YLYTLT3 \
           --sample=CMU_CV_B3-1_GEX \
           --create-bam=true \
           --localcores=64 \
           --localmem=256 \
           --tenx-cloud-token-path=/home/warunyoo.p/.config/txg/credentials

cellranger count --id=visit1-2 \
           --transcriptome=/common/db/human_ref/cellranger/refdata-gex-GRCh38-2024-A \
           --fastqs=/home/warunyoo.p/o250031_COVID68/HN00254955/HN00254955_10X_RawData_Outs/CMU_CV_B3-2_GEX/22YLYTLT3 \
           --sample=CMU_CV_B3-2_GEX \
           --create-bam=true \
           --localcores=64 \
           --localmem=256 \
           --tenx-cloud-token-path=/home/warunyoo.p/.config/txg/credentials

cellranger count --id=visit2-1 \
           --transcriptome=/common/db/human_ref/cellranger/refdata-gex-GRCh38-2024-A \
           --fastqs=/home/warunyoo.p/o250031_COVID68/HN00254955/HN00254955_10X_RawData_Outs/CMU_CV_B4-1_GEX/22YLYTLT3 \
           --sample=CMU_CV_B4-1_GEX \
           --create-bam=true \
           --localcores=64 \
           --localmem=256 \
           --tenx-cloud-token-path=/home/warunyoo.p/.config/txg/credentials

cellranger count --id=visit2-2 \
           --transcriptome=/common/db/human_ref/cellranger/refdata-gex-GRCh38-2024-A \
           --fastqs=/home/warunyoo.p/o250031_COVID68/HN00254955/HN00254955_10X_RawData_Outs/CMU_CV_B4-2_GEX/22YLYTLT3 \
           --sample=CMU_CV_B4-2_GEX \
           --create-bam=true \
           --localcores=64 \
           --localmem=256 \
           --tenx-cloud-token-path=/home/warunyoo.p/.config/txg/credentials