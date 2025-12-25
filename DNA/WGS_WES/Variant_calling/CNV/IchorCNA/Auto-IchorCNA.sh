#!/bin/bash
#SBATCH --account=o250000               
#SBATCH --job-name=IchorCNA          
#SBATCH --partition=compute             
#SBATCH --time=05:00:00               
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G               
#SBATCH --ntasks-per-node=1            
#SBATCH --output=Ichor_%j_output.log          
#SBATCH --error=Ichor_%j_error.log            

ml purge
ml mamba/1.5.8
ml ichorCNA/0.2.0

GC_19='/apps/ichorCNA/0.2.0/inst/extdata/gc_hg19_1000kb.wig' 
MAP_19='/apps/ichorCNA/0.2.0/inst/extdata/map_hg19_1000kb.wig'
GC_38='/apps/ichorCNA/0.2.0/inst/extdata/gc_hg38_1000kb.wig' 
MAP_38='/apps/ichorCNA/0.2.0/inst/extdata/map_hg38_1000kb.wig'
GC_19_WES='/apps/ichorCNA/0.2.0/inst/extdata/gc_hg19_1000kb.wig' 
MAP_19_WES='/apps/ichorCNA/0.2.0/inst/extdata/map_hg19_1000kb.wig'
GC_38_WES='/apps/ichorCNA/0.2.0/inst/extdata/gc_hg38_10kb.wig' 
MAP_38_WES='/apps/ichorCNA/0.2.0/inst/extdata/map_hg38_10kb.wig'
EXONs_38='/apps/ichorCNA/0.2.0/inst/extdata/Exon_regions_hg38.v1.bed'

IchorCNA_script='/apps/ichorCNA/0.2.0/scripts/runIchorCNA.R'

chromosome='chr'

if [[ "${chromosome}" = "num" ]];then
    chrom='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X'
elif [[ "${chromosome}" = "chr" ]];then
    chrom='chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX'
fi

for i in `cat /project/o250003_CFBile/WGS_seq/sample_all.txt`
do
    readCounter --window 1000000 \
        --quality 20 \
        --chromosome ${chrom} \
        /project/o250003_CFBile/WGS_seq/Preprocess/${i}/${i}_recal.bam > /project/o250003_CFBile/WGS_seq/Preprocess/${i}/${i}.wig 

    IchorCNA_WGS_hg19() {
        Rscript ${IchorCNA_script} --id ${i} \
        --WIG /project/o250003_CFBile/WGS_seq/Preprocess/${i}/${i}.wig --ploidy "c(2)" \
        --normal "c(0.5,0.6,0.7,0.8,0.9)" \
        --gcWig ${GC_19} \
        --mapWig ${MAP_19} \
        --centromere /apps/ichorCNA/0.2.0/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
        --includeHOMD False --chrs "c(1:22)" --chrTrain "c(1:22)" --genomeBuild "hg19" --genomeStyle "UCSC" \
        --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
        --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir /project/o250003_CFBile/WGS_seq/Preprocess/${i}
    }
    IchorCNA_WGS_hg38() {
        Rscript ${IchorCNA_script} --id ${i} \
        --WIG /project/o250003_CFBile/WGS_seq/Preprocess/${i}/${i}.wig --ploidy "c(2)" \
        --normal "c(0.5,0.6,0.7,0.8,0.9)" \
        --gcWig ${GC_38} \
        --mapWig ${MAP_38} \
        --centromere /apps/ichorCNA/0.2.0/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
        --includeHOMD False --chrs "c(1:22)" --chrTrain "c(1:22)" --genomeBuild "hg38" --genomeStyle "UCSC" \
        --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
        --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir /project/o250003_CFBile/WGS_seq/Preprocess/${i}
    }

    # IchorCNA_WGS_hg19
    IchorCNA_WGS_hg38

done


