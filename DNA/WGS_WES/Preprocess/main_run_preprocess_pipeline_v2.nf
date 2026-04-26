#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// =============================================================================
// Pipeline   : HNC QC Preprocessing Pipeline (Nextflow DSL2)
// Description: FastQC, Fastp, Parabricks fq2bam/applybqsr, Samtools, MosDepth, MultiQC
// Converted from: run_preprocess_pipeline_HN00270111.sh
// =============================================================================

// --- Params ---
params.samples_sheet   = 'NO_FILE'  // CSV with columns: ID, R1, R2
params.out_dir         = 'NO_DIR'   // Base output directory
params.fasta           = System.getenv('REF_FA') ?: '/common/db/human_ref/hg38/parabricks/Homo_sapiens_assembly38.fasta'
params.known_sites_dir = '/common/db/human_ref/hg38/v0'
params.exome_bed       = '/project/o260003_CRU_BI/human_ref/exome_calling_regions.v1.interval_list.bed'
params.qc_script       = 'NO_FILE'  // Path to qc_script_v1.2.py
params.sample_regex    = '^([^_]+(?:_[^_]+){5})'

// Container paths (override via env or CLI)
params.pb_container       = System.getenv('PB_SIF')       ?: '/common/sif/clara-parabricks/4.4.0.sif'
params.mosdepth_container = System.getenv('MOSDEPTH_SIF') ?: '/common/sif/mosdepth/Mosdepth_0.3.10--h4e814b3_1.sif'

// GPU settings
params.num_gpus        = 1
params.slurm_account   = System.getenv('SLURM_ACCOUNT') ?: 'o250022'

// =============================================================================
// PROCESS 1: FASTQC_RAW
// =============================================================================
process FASTQC_RAW {
    tag "${sample_id}"
    label 'cpu_low'
    publishDir "${params.out_dir}/reports/fastqc_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "*_fastqc.html", emit: html
    path "*_fastqc.zip",  emit: zip

    script:
    """
    module load fastqc
    fastqc -t ${task.cpus} -q ${r1} ${r2}
    """
}

// =============================================================================
// PROCESS 2: FASTP
// =============================================================================
process FASTP {
    tag "${sample_id}"
    label 'cpu_med'
    publishDir "${params.out_dir}/reports", mode: 'copy', pattern: "*.{json,html}"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_R1_clean.fastq.gz"), path("${sample_id}_R2_clean.fastq.gz"), emit: reads
    path("${sample_id}.json"), emit: json_report
    path("${sample_id}.html"), emit: html_report

    script:
    """
    module load fastp
    fastp \
        --thread ${task.cpus} \
        --in1 ${r1} \
        --in2 ${r2} \
        --out1 ${sample_id}_R1_clean.fastq.gz \
        --out2 ${sample_id}_R2_clean.fastq.gz \
        --json ${sample_id}.json \
        --html ${sample_id}.html \
        --detect_adapter_for_pe \
        --trim_poly_g
    """
}

// =============================================================================
// PROCESS 3: FASTQC_CLEAN
// =============================================================================
process FASTQC_CLEAN {
    tag "${sample_id}"
    label 'cpu_low'
    publishDir "${params.out_dir}/reports/fastqc_clean", mode: 'copy'

    input:
    tuple val(sample_id), path(r1_clean), path(r2_clean)

    output:
    path "*_fastqc.html", emit: html
    path "*_fastqc.zip",  emit: zip

    script:
    """
    module load fastqc
    fastqc -t ${task.cpus} -q ${r1_clean} ${r2_clean}
    """
}

// =============================================================================
// PROCESS 4: FQ2BAM (GPU — Parabricks via Apptainer)
// =============================================================================
process FQ2BAM {
    tag "${sample_id}"
    label 'gpu'
    container params.pb_container

    input:
    tuple val(sample_id), path(r1_clean), path(r2_clean)
    val(interval)

    output:
    tuple val(sample_id),
          path("${sample_id}_dedup_sorted.bam"),
          path("${sample_id}_dedup_sorted.bam.bai"),
          path("${sample_id}-recal.txt"), emit: bam_recal
    path("${sample_id}-qc-metrics/*"), emit: qc_metrics

    script:
    """
    pbrun fq2bam \
        --num-gpus ${params.num_gpus} \
        --ref ${params.fasta} \
        --in-fq ${r1_clean} ${r2_clean} \
        --knownSites ${params.known_sites_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        --knownSites ${params.known_sites_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz \
        --knownSites ${params.known_sites_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
        --interval ${params.known_sites_dir}/${interval} \
        --read-group-sm ${sample_id} \
        --read-group-lb ${sample_id} \
        --read-group-pl Illumina \
        --read-group-id-prefix ${sample_id} \
        --out-bam ${sample_id}_dedup_sorted.bam \
        --out-qc-metrics-dir ${sample_id}-qc-metrics/ \
        --out-recal-file ${sample_id}-recal.txt \
        --tmp-dir .
    """
}

// =============================================================================
// PROCESS 5: APPLYBQSR (GPU — Parabricks via Apptainer)
// =============================================================================
process APPLYBQSR {
    tag "${sample_id}"
    label 'gpu'
    container params.pb_container

    input:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai), path(recal_txt)
    val(interval)

    output:
    tuple val(sample_id), path("${sample_id}_recal.bam"), path("${sample_id}_recal.bam.bai"), emit: recal_bam
    tuple val(sample_id), path(dedup_bam), path(dedup_bai), emit: dedup_bam

    script:
    """
    pbrun applybqsr \
        --num-gpus ${params.num_gpus} \
        --ref ${params.fasta} \
        --in-bam ${dedup_bam} \
        --in-recal-file ${recal_txt} \
        --interval ${params.known_sites_dir}/${interval} \
        --out-bam ${sample_id}_recal.bam \
        --tmp-dir .
    """
}

// =============================================================================
// PROCESS 6: BAM_POSTPROCESS (flagstat + BAM-to-CRAM + cleanup)
// =============================================================================
process BAM_POSTPROCESS {
    tag "${sample_id}"
    label 'cpu_med'
    publishDir "${params.out_dir}/02_PREPROCESS/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(recal_bam), path(recal_bai)
    tuple val(sample_id2), path(dedup_bam), path(dedup_bai)

    output:
    path("${sample_id}_flagstat.txt"), emit: flagstat
    tuple path("${sample_id}_dedup_sorted.cram"), path("${sample_id}_dedup_sorted.cram.crai"), emit: cram

    script:
    """
    module load samtools

    samtools flagstats -@ ${task.cpus} ${recal_bam} > ${sample_id}_flagstat.txt

    samtools view -@ ${task.cpus} -C -T ${params.fasta} \
        -o ${sample_id}_dedup_sorted.cram ${dedup_bam}

    samtools index -@ ${task.cpus} -b ${sample_id}_dedup_sorted.cram
    """
}

// =============================================================================
// PROCESS 7: MOSDEPTH (container, CPU)
// =============================================================================
process MOSDEPTH {
    tag "${sample_id}"
    label 'cpu_med'
    container params.mosdepth_container
    publishDir "${params.out_dir}/02_PREPROCESS/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(recal_bam), path(recal_bai)
    val(seq_platform)

    output:
    tuple val(sample_id), path("${sample_id}.regions.bed.gz"), emit: regions
    path("${sample_id}.mosdepth.summary.txt"), emit: summary

    script:
    def by_param = seq_platform == 'WGS' ? '--by 1000' : "--by ${params.exome_bed}"
    """
    mosdepth -n --fast-mode -t ${task.cpus} \
        ${by_param} \
        ${sample_id} ${recal_bam}
    """
}

// =============================================================================
// PROCESS 8: CALC_REGION_DEPTH
// =============================================================================
process CALC_REGION_DEPTH {
    tag "${sample_id}"
    label 'cpu_low'

    input:
    tuple val(sample_id), path(regions_bed_gz)

    output:
    path("${sample_id}_depth.tsv"), emit: depth_line

    script:
    """
    zcat ${regions_bed_gz} \
      | awk '\$1 ~ /^chr[0-9]+\$/ { n=substr(\$1,4)+0; if(n>=1 && n<=22) print \$5 }' \
      | sort -g \
      | awk '{sum+=\$1; vals[NR]=\$1}
        END {
          if(NR==0){printf "${sample_id}\\tNA\\tNA\\n"; exit}
          mean=sum/NR
          if(NR%2==1) median=vals[int(NR/2)+1]
          else median=(vals[NR/2]+vals[NR/2+1])/2
          printf "${sample_id}\\t%.4f\\t%.4f\\n", mean, median
        }' > ${sample_id}_depth.tsv
    """
}

// =============================================================================
// PROCESS 9: MULTIQC_AND_REPORT
// =============================================================================
process MULTIQC_AND_REPORT {
    label 'cpu_low'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path('reports/*')
    path(region_depth_tsv)

    output:
    path("multiqc_report.html"), emit: multiqc_html
    path("multiqc_data"), emit: multiqc_data
    path("QC_SUMMARY.xlsx"), emit: qc_summary

    script:
    """
    module load multiqc
    multiqc reports/ -o .

    python ${params.qc_script} \
        -i multiqc_data/multiqc_data.json \
        -o QC_SUMMARY.xlsx \
        -p '${params.sample_regex}' \
        -d ${region_depth_tsv}
    """
}

// =============================================================================
// HELPER: Determine sequencing platform from sample name
// =============================================================================
def get_seq_platform(sample_id) {
    if (sample_id.contains('WGS')) {
        return 'WGS'
    } else if (sample_id.contains('WES')) {
        return 'WES'
    } else {
        error "Cannot determine platform (WGS/WES) for sample: ${sample_id}"
    }
}

def get_interval(sample_id) {
    if (sample_id.contains('WGS')) {
        return 'wgs_calling_regions.hg38.interval_list'
    } else if (sample_id.contains('WES')) {
        return 'exome_calling_regions.v1.interval_list'
    } else {
        error "Cannot determine interval for sample: ${sample_id}"
    }
}

// =============================================================================
// WORKFLOW
// =============================================================================
workflow {
    // Validate required params
    if (params.samples_sheet == 'NO_FILE') { error "Please provide --samples_sheet" }
    if (params.out_dir == 'NO_DIR')        { error "Please provide --out_dir" }
    if (params.qc_script == 'NO_FILE')     { error "Please provide --qc_script" }

    if (!file(params.samples_sheet).exists())
        error "Samplesheet not found: ${params.samples_sheet}"

    // Create sample channel from CSV: ID, R1, R2
    samples_ch = Channel
        .fromPath(params.samples_sheet, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def sample = row.ID?.toString()?.trim()
            def r1     = row.R1?.toString()?.trim()
            def r2     = row.R2?.toString()?.trim()

            if (!sample || !r1 || !r2)
                error "Samplesheet must be CSV with non-empty columns: ID, R1, R2"

            tuple(sample, file(r1, checkIfExists: true), file(r2, checkIfExists: true))
        }

    // 1. FastQC on raw reads
    FASTQC_RAW(samples_ch)

    // 2. Fastp trimming
    FASTP(samples_ch)

    // 3. FastQC on cleaned reads
    FASTQC_CLEAN(FASTP.out.reads)

    // 4. Parabricks fq2bam — need interval per sample
    fq2bam_input = FASTP.out.reads.map { sample_id, r1, r2 ->
        tuple(sample_id, r1, r2)
    }
    fq2bam_interval = FASTP.out.reads.map { sample_id, r1, r2 ->
        get_interval(sample_id)
    }
    FQ2BAM(fq2bam_input, fq2bam_interval)

    // 5. ApplyBQSR
    applybqsr_interval = FQ2BAM.out.bam_recal.map { sample_id, bam, bai, recal ->
        get_interval(sample_id)
    }
    APPLYBQSR(FQ2BAM.out.bam_recal, applybqsr_interval)

    // 6. BAM post-processing (flagstat, BAM-to-CRAM, cleanup)
    BAM_POSTPROCESS(APPLYBQSR.out.recal_bam, APPLYBQSR.out.dedup_bam)

    // 7. MosDepth coverage
    mosdepth_platform = APPLYBQSR.out.recal_bam.map { sample_id, bam, bai ->
        get_seq_platform(sample_id)
    }
    MOSDEPTH(APPLYBQSR.out.recal_bam, mosdepth_platform)

    // 8. Calculate region depth per sample
    CALC_REGION_DEPTH(MOSDEPTH.out.regions)

    // Collect all depth lines into a single TSV with header
    region_depth_tsv = CALC_REGION_DEPTH.out.depth_line
        .collectFile(
            name: 'mosdepth_region_depth.tsv',
            storeDir: "${params.out_dir}/reports",
            seed: 'Sample\tMean_Region_Depth_chr1_22\tMedian_Region_Depth_chr1_22\n',
            newLine: false
        )

    // 9. MultiQC + QC report — collect all reports
    all_reports = FASTQC_RAW.out.html.mix(FASTQC_RAW.out.zip)
        .mix(
            FASTP.out.json_report,
            FASTP.out.html_report,
            FASTQC_CLEAN.out.html,
            FASTQC_CLEAN.out.zip,
            FQ2BAM.out.qc_metrics,
            BAM_POSTPROCESS.out.flagstat,
            MOSDEPTH.out.summary
        )
        .collect()

    MULTIQC_AND_REPORT(all_reports, region_depth_tsv)
}
