#!/bin/bash

#SBATCH --account=o250022               
#SBATCH --job-name=fq2bam  
#SBATCH --partition=gpu
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=20
#SBATCH --gpus=1
#SBATCH --mem=200G
#SBATCH --output=fq2bam_output_%j.log          
#SBATCH --error=fq2bam_error_%j.log


inpdir='/project/o250003_CFBile/Target_Seq/02_Macrogen_Prep/01_CLEAN_QC'
outdir='/project/o250003_CFBile/Target_Seq/02_Macrogen_Prep/02_PREPROCESS'

for ID in `cat /project/o250003_CFBile/Target_Seq/02_Macrogen_Prep/00_RAW_FASTQ/Round1_id.txt`
do
    # Load required modules
    ml purge
    ml fastqc/0.12.1
    ml apptainer
    ml samtools

    # gatkdir=${outdir}/${ID}/04_BQSR
    noutdir=${outdir}/${ID}/01_FQ2BAM

    mkdir ${noutdir}

    # FQ2BAM
    apptainer exec --nv \
        -B /common/db/human_ref/hg38/parabricks:/Ref_hg38 \
        -B /common/db/human_ref/hg38/v0:/Ref_hg38_v0 \
        -B ${inpdir}:/FASTQdir \
        -B ${noutdir}:/Workdir \
        /common/sif/clara-parabricks/4.4.0.sif pbrun \
            fq2bam \
            --ref /Ref_hg38/Homo_sapiens_assembly38.fasta \
            --in-fq /FASTQdir/${ID}_1_cleaned.fastq.gz /FASTQdir/${ID}_2_cleaned.fastq.gz \
            --knownSites /Ref_hg38_v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
            --knownSites /Ref_hg38_v0/Homo_sapiens_assembly38.known_indels.vcf.gz \
            --knownSites /Ref_hg38_v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
            --read-group-id-prefix ${ID} \
            --out-bam /Workdir/${ID}_dedup_sorted.bam \
            --out-qc-metrics-dir /Workdir/${ID}-qc-metrics/ \
            --out-recal-file /Workdir/${ID}-recal.txt \
            --tmp-dir /Workdir

    # ApplyBQSR
    apptainer exec --nv \
        -B /common/db/human_ref/hg38/parabricks:/Ref_hg38 \
        -B ${noutdir}:/InputDir \
        -B ${noutdir}:/OutputDir \
        -B /common/db/human_ref/hg38/v0:/Ref_hg38_v0 \
        /common/sif/clara-parabricks/4.4.0.sif pbrun \
            applybqsr \
            --ref /Ref_hg38/Homo_sapiens_assembly38.fasta \
            --in-bam /InputDir/${ID}_dedup_sorted.bam \
            --in-recal-file /InputDir/${ID}-recal.txt  \
            --out-bam /OutputDir/${ID}_recal.bam \
            --tmp-dir /OutputDir

    samtools index -@ 20 -b ${noutdir}/${ID}_recal.bam

    # MosDepth For WGS
    apptainer run --nv \
        -B ${noutdir}:/InputDir \
        -B ${noutdir}:/OutputDir \
        -B /project/o250003_CFBile/Target_Seq/02_Macrogen_Prep/src:/Onco_pan \
        /home/songphon.sut/SIF/Mosdepth_0.3.10--h4e814b3_1.sif \
        mosdepth -n --fast-mode -t 20 \
        --by 1000 \
        /OutputDir/${ID}_WGS_depth /InputDir/${ID}_recal.bam

    # MosDepth For WES or Target
    # apptainer run --nv \
    #     -B ${noutdir}:/InputDir \
    #     -B ${noutdir}:/OutputDir \
    #     -B /project/o250003_CFBile/Target_Seq/02_Macrogen_Prep/src:/Onco_pan \
    #     /home/songphon.sut/SIF/Mosdepth_0.3.10--h4e814b3_1.sif \
    #     mosdepth -n --fast-mode -t 20 \
    #     --by /Onco_pan/onco_680_DNA_hg38.sorted_merged.bed \
    #     /OutputDir/${ID}_WGS_depth /InputDir/${ID}_recal.bam

    # Flagstat
    samtools flagstats -@ 20 ${noutdir}/${ID}_recal.bam > ${noutdir}/${ID}_flagstat.txt

done

