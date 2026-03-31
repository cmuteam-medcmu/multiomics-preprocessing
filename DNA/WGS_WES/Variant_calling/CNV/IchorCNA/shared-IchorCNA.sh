#!/bin/bash
#SBATCH --account=o260002              
#SBATCH --job-name=IchorCNA           
#SBATCH --partition=compute              
#SBATCH --time=05:00:00                
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G                
#SBATCH --ntasks-per-node=1             
#SBATCH --output=IchorCNA_%j.out         
#SBATCH --error=IchorCNA_%j.err          

set -euo pipefail

# 1. DYNAMIC PARAMETERS
readonly BASE_DIR="/project/o260002_CFSTS69/Tissue/WES-500x/20260204_HN00265200_EXO_Report"
readonly SAMPLES_SHEET="${BASE_DIR}/sample_sheet.csv"

# 2. STATIC PARAMETERS
readonly GENOME_BUILD="hg38"
readonly ICHOR_DATA="/apps/ichorCNA/0.2.0/inst/extdata"
readonly ICHOR_SCRIPT="/apps/ichorCNA/0.2.0/scripts/runIchorCNA.R"
readonly CHROM_LIST='chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX'


# 3. FUNCTIONS
run_readcounter() {
    local sample_id=$1 local bam_input=$2 local output_wig=$3
    ml purge; ml mamba/1.5.8 ichorCNA/0.2.0
    echo "LOG: Running readCounter for ${sample_id}..."
    readCounter --window 1000000 --quality 20 --chromosome "${CHROM_LIST}" "${bam_input}" > "${output_wig}"
}

run_ichorcna() {
    local sample_id=$1 local wig_input=$2 local output_dir=$3 local build=$4

    echo "LOG: [$(date)] Running IchorCNA (${build}) for ${sample_id}..."
    ml purge; ml mamba/1.5.8 ichorCNA/0.2.0

    local gc_wig="${ICHOR_DATA}/gc_${build}_1000kb.wig"
    local map_wig="${ICHOR_DATA}/map_${build}_1000kb.wig"
    local centromere="${ICHOR_DATA}/GRCh37.p13_centromere_UCSC-gapTable.txt"
    [[ "$build" == "hg38" ]] && centromere="${ICHOR_DATA}/GRCh38.GCA_000001405.2_centromere_acen.txt"

    Rscript "${ICHOR_SCRIPT}" --id "${sample_id}" --WIG "${wig_input}" \
        --normal "c(0.5,0.6,0.7,0.8,0.9)" --gcWig "${gc_wig}" --mapWig "${map_wig}" \
        --centromere "${centromere}" --includeHOMD False --chrs "c(1:22)" --chrTrain "c(1:22)" \
        --genomeBuild "${build}" --genomeStyle "UCSC" --estimateNormal True \
        --estimatePloidy True --estimateScPrevalence True --scStates "c(1,3)" \
        --txnE 0.9999 --txnStrength 10000 --outDir "${output_dir}"
}


# 4. MAIN LOGIC
main() {
    echo "LOG: Starting IchorCNA Pipeline"
    
    {
        read -r header # Skip CSV header (ID,R1,R2)
        
        while IFS=',' read -r sample_id r1 r2 || [[ -n "$sample_id" ]]; do
            [[ -z "$sample_id" ]] && continue

            local sample_outdir="${BASE_DIR}/04_CNV/${sample_id}"
            local bam_input="${BASE_DIR}/02_PREPROCESS/${sample_id}/${sample_id}_recal.bam"
            local wig_file="${sample_outdir}/${sample_id}.wig"

            # Proactive check: If .bam not found, try .cram from previous step
            [[ ! -f "$bam_input" ]] && bam_input="${BASE_DIR}/02_PREPROCESS/${sample_id}/${sample_id}_recal.cram"

            mkdir -p "${sample_outdir}"

            if [[ -f "$bam_input" ]]; then
                run_readcounter "$sample_id" "$bam_input" "$wig_file"
                run_ichorcna "$sample_id" "$wig_file" "$sample_outdir" "$GENOME_BUILD"
            else
                echo "SKIP: File not found for ${sample_id} at ${bam_input}"
            fi
        done
    } < "$SAMPLES_SHEET"

    echo "SUCCESS: Pipeline completed."
}

main "$@"