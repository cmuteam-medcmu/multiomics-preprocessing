# Somatic Variant Calling — Clara Parabricks / Mutect2

GPU-accelerated somatic variant calling using [NVIDIA Clara Parabricks](https://docs.nvidia.com/clara/parabricks/) on a SLURM HPC cluster. This module runs paired tumor–normal variant calling with GATK Mutect2 via two distinct entry points depending on available input data.

---

## Overview

Two pipeline modes are supported:

| Script | Input | Entry Point | Use Case |
|---|---|---|---|
| `clara_pbrun_somatic_WES-WES.sh` | Tumor (WES) FASTQs + Normal (WES) FASTQs | `pbrun somatic` | Full alignment + BQSR + calling in one step |
| `clara_pbrun_somatic_WGS-WES.sh` | Tumor (WES) FASTQs + Normal (WGS) FASTQs | `pbrun somatic` | Mixed-modality: WGS Normal with exome Tumor |
| `clara_pbrun_mutectcaller_WES.sh` | Pre-aligned Tumor BAM + Normal BAM | `pbrun mutectcaller` | Variant calling only from existing recalibrated BAMs |

> **When to use `somatic` vs `mutectcaller`:**
> Use the `somatic` scripts when starting from FASTQs — they run BWA-MEM2 alignment, BQSR, and Mutect2 in a single optimized GPU pass. Use `mutectcaller` when you already have recalibrated BAMs (e.g., output from a previous `somatic` run or another preprocessing pipeline).

---

## Pipeline Diagram

```
FASTQs (Tumor + Normal)
        │
        ▼
  [pbrun somatic]                          Pre-aligned BAMs
  ┌─────────────────────────────┐                │
  │  BWA-MEM2 Alignment         │                ▼
  │  Base Quality Score Recal.  │    [pbrun mutectcaller]
  │  Mutect2 Variant Calling    │    Mutect2 Variant Calling
  └─────────────┬───────────────┘                │
                │                                │
                ▼                                ▼
        Recalibrated BAMs             Unfiltered VCF (.vcf.gz)
        Recal Tables (.txt)
        Unfiltered VCF (.vcf.gz)
```

---

## Requirements

- **SLURM** workload manager
- **Apptainer** (Singularity-compatible)
- **Clara Parabricks** SIF image: `/common/sif/clara-parabricks/4.5.1-1.sif`
- **NVIDIA GPU** (1× GPU per job)
- **bcftools** (for post-run variant counting)
- Reference genome: **hg38** (GRCh38)

---

## Reference Files Required

| File | Description |
|---|---|
| `Homo_sapiens_assembly38.fasta` | hg38 reference genome (Parabricks-indexed) |
| `1000G_phase1.snps.high_confidence.hg38.vcf.gz` | Known SNP sites for BQSR |
| `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` | Known indel sites for BQSR |
| Exon capture BED file | Target regions BED for your WES kit |

---

## Directory Structure

Before running, configure the following path variables inside each script:

```bash
PROJECT_DIR=/path_to_project_DIR
REF_DIR=/path_to_reference_DIR           # hg38 Parabricks reference
KNOWN_SITES_DIR=/path_to_KNOWN_SITES_DIR # BQSR known sites VCFs
FASTQ_DIR=${PROJECT_DIR}/path_to_fp.fastq_DIR
BAM_DIR=${PROJECT_DIR}/path_to_BAM_storage
BED=/path_to_Exon_BED_file.bed
```

Expected output layout (auto-created by the scripts):

```
Somatic_callers/
└── <PATIENT_ID>/
    ├── vcf/
    │   └── mutect2/
    │       └── <PATIENT_ID>.unfiltered.vcf.gz
    ├── recal_tables/
    │   ├── <PATIENT_ID>_tumor_WES.recal.txt
    │   └── <PATIENT_ID>_normal_W[E/G]S.recal.txt
    └── tmp/                  # removed automatically on success
```

---

## Usage

All scripts accept a single positional argument: the patient/sample ID.

```bash
# Full pipeline from FASTQs — WES Tumor vs WES Normal
sbatch clara_pbrun_somatic_WES-WES.sh <PATIENT_ID>

# Full pipeline from FASTQs — WES Tumor vs WGS Normal
sbatch clara_pbrun_somatic_WGS-WES.sh <PATIENT_ID>

# Calling only from pre-aligned BAMs
sbatch clara_pbrun_mutectcaller_WES.sh <PATIENT_ID>
```

---

## FASTQ Naming Convention

The `somatic` scripts distinguish tumor from normal FASTQs by filename suffix. Update the patterns in the script to match your data:

```bash
# Tumor (WES) — suffix _T_
FASTQ_T1=${FASTQ_DIR}/${PATIENT}_T_read_1.fp.fastq.gz
FASTQ_T2=${FASTQ_DIR}/${PATIENT}_T_read_2.fp.fastq.gz

# Normal (WES or WGS) — suffix _N_
FASTQ_N1=${FASTQ_DIR}/${PATIENT}_N_read_1.fp.fastq.gz
FASTQ_N2=${FASTQ_DIR}/${PATIENT}_N_read_2.fp.fastq.gz
```

> ⚠️ **Important:** In the original scripts both tumor and normal FASTQ variables pointed to identical paths. This has been corrected — verify that your actual file names differ between tumor and normal samples.

---

## Output

| File | Description |
|---|---|
| `<PATIENT>.unfiltered.vcf.gz` | Raw Mutect2 calls (not yet filtered) |
| `<PATIENT>_tumor_WES.recal.bam` | Recalibrated tumor BAM (`somatic` mode only) |
| `<PATIENT>_normal_W[E/G]S.recal.bam` | Recalibrated normal BAM (`somatic` mode only) |
| `*.recal.txt` | BQSR recalibration tables (`somatic` mode only) |

The unfiltered VCF is the input for downstream filtering steps (e.g., `FilterMutectCalls`, PON filtering, annotation).

---

## Notes and Known Limitations

- **Unfiltered output only.** These scripts produce raw Mutect2 calls. Downstream filtering (e.g., `pbrun filtermutect` or GATK `FilterMutectCalls`) is a separate step not included here.
- **WGS Normal with exome intervals.** In the `WGS-WES` mode, the WGS Normal FASTQ is aligned and the `--interval` BED is applied at the calling step. Coverage outside the BED is discarded. This is an appropriate approach for hybrid protocols.
- **`mutectcaller` vs `somatic` resource allocation.** The `mutectcaller` script requests more CPUs (47) and memory (240 GB) than the `somatic` scripts (24 CPUs, 120 GB) because it performs no alignment but runs a CPU-intensive Mutect2 pass. Adjust `--cpus-per-task` and `--mem` to match your cluster's available nodes.
- **Temp directory cleanup.** On successful completion, `${TMP_DIR}` is removed automatically. If the job fails, the temp directory is preserved for debugging.

---

## SLURM Resource Summary

| Script | CPUs | Memory | GPU | Time Limit |
|---|---|---|---|---|
| `clara_pbrun_somatic_WES-WES.sh` | 24 | 120 GB | 1 | 24 h |
| `clara_pbrun_somatic_WGS-WES.sh` | 24 | 120 GB | 1 | 24 h |
| `clara_pbrun_mutectcaller_WES.sh` | 47 | 240 GB | 1 | 24 h |
