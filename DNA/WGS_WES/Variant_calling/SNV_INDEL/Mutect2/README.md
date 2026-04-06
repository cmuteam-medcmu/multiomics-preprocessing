# 🧬 Somatic Variant Calling Pipeline (Mutect2-based)

## Overview

This pipeline performs **tumor–normal somatic variant calling** using **GATK Mutect2**, followed by standard GATK filtering and optional **custom post-filtering** for improved variant confidence.

This pipeline assumes **paired tumor–normal samples** for somatic variant detection.

It is designed as a **shared, reproducible workflow** that:

* supports multiple samples via a sample sheet
* separates CPU and GPU steps
* allows flexible filtering strategies
* supports both **parallel** and **serial (one-case-at-a-time)** execution

---

## Pipeline Structure

The pipeline consists of three main steps:

### 1. GATK Preprocessing Tables

* Tool: `GetPileupSummaries`, `CalculateContamination`
* Input: tumor + normal BAM
* Output:

  * pileup summaries
  * contamination table
  * segmentation table

---

### 2. Somatic Variant Calling

* Tool: `Mutect2`

  * GPU-accelerated version via Parabricks (recommended)
  * or standard GATK CPU implementation
* Input:

  * tumor BAM
  * normal BAM
  * contamination estimates
* Output:

  * unfiltered somatic VCF

---

### 3. Post-filtering

* Tool: `FilterMutectCalls` + optional custom filtering
* Output:

  * GATK-filtered VCF
  * optional stricter high-confidence VCF

---

## Sample Sheet Format

The pipeline uses a CSV sample sheet:

```csv
case_id,tumor_sample_id,normal_sample_id,tumor_bam,normal_bam
LK00019,20260324_WES500X_ILL_LK00019_BMA_Dx,20260324_WES100X_ILL_LK00019_BC_Dx,/path/tumor.bam,/path/normal.bam
```

### Columns

| Column           | Description                                  |
| ---------------- | -------------------------------------------- |
| case_id          | Unique patient/case ID                       |
| tumor_sample_id  | Tumor sample identifier (used in read group) |
| normal_sample_id | Normal sample identifier                     |
| tumor_bam        | Path to tumor BAM                            |
| normal_bam       | Path to normal BAM                           |

---

## How to Run

```bash
bash submit_mutect2_pipeline.sh
```

---

## Execution Model

### Within a case

* Step 1 and Step 2 run in parallel
* Step 3 waits for both

### Across cases (controlled by `SERIAL_BY_CASE`)

#### Serial mode

```
Case1 → Case2 → Case3
```

#### Parallel mode

```
Case1
Case2
Case3   (independent submission)
```

---

## Configuration

### Key parameters (in submit script)

```bash
readonly FILTER_MODE="all"
readonly SERIAL_BY_CASE="yes"
readonly SLURM_ACCOUNT="o260003"
```

---

### FILTER_MODE options

| Mode        | Description                                                 |
| ----------- | ----------------------------------------------------------- |
| `gatk_pass` | Output variants with `FILTER=PASS` from GATK only           |
| `custom_qc` | Apply custom QC filters and output high-confidence variants |
| `all`       | Output both GATK PASS and custom QC-filtered results        |

---

### SERIAL_BY_CASE

| Value | Behavior                            |
| ----- | ----------------------------------- |
| `yes` | Run one case at a time (sequential) |
| `no`  | Submit all cases in parallel        |

---

## Filtering Strategy

### 1. GATK Filtering

Variants are first filtered using **GATK FilterMutectCalls**, which assigns:

```
FILTER = PASS or FAIL (e.g., strand_bias, contamination)
```

---

### 2. Custom QC Filtering (optional)

Additional filters are applied using `bcftools`:

```
Tumor DP ≥ 20
AND Normal DP ≥ 10
AND Tumor alt AD ≥ 5
```

These thresholds are chosen to ensure:

* sufficient tumor coverage for reliable variant detection
* adequate alternate allele support to reduce sequencing noise
* sufficient normal coverage to support somatic classification

---

### Important Logic

* Variants that fail any condition are **flagged**:

  * `LOW_T_DP`
  * `LOW_N_DP`
  * `LOW_T_AD`

* The pipeline then keeps only variants where:

```
FILTER == PASS
```

### Key Interpretation

This does **not redefine PASS**. Instead, it retains only variants that:

```
GATK PASS
AND do not trigger any custom filter flags
```

Therefore, the final high-confidence variants satisfy:

```
GATK PASS
AND Tumor DP ≥ 20
AND Normal DP ≥ 10
AND Tumor alt AD ≥ 5
```

---

## Target Regions (WES)

Target intervals (e.g., BED file) are defined within the pipeline scripts.

Ensure that:

* the interval file matches the capture kit used (e.g., Agilent SureSelect V8)
* tumor and normal BAMs are aligned consistently to the same reference and intervals

---

## Output Structure

```
output/
├── LK00019/
│   ├── gatk_table/
│   ├── mutect2/
│   ├── postfilter/
│   ├── vcf/
│   │   ├── *.filtered.vcf.gz              # GATK output
│   │   ├── *.filtered.PASS.vcf.gz         # GATK PASS only
│   │   ├── *.filtered.PASS.QUALITY.vcf.gz # GATK PASS + custom QC
```

---

## Logging

Logs are organized per case:

```
logs/
├── LK00019/
│   ├── gatk_table/
│   ├── mutect2/
│   ├── postfilter/
```

Each job produces:

* `*.log`
* `*.err`

---

## Resource Usage

Each step is executed as an independent Slurm job.

### Example configuration

| Step        | CPU       | Memory | GPU |
| ----------- | --------- | ------ | --- |
| GATK tables | 4 cores   | 32 GB  | ❌   |
| Mutect2     | 40+ cores | 150 GB | ✅   |
| Postfilter  | 2 cores   | 8 GB   | ❌   |

These resources are **per sample**, not per sample sheet.

---

## Design Principles

* **Reproducibility**: standardized sample sheet input
* **Modularity**: each step is an independent script
* **Scalability**: supports both HPC parallel and sequential modes
* **Flexibility**: configurable filtering strategy
* **Clarity**: explicit separation of GATK vs custom QC filtering

---

## Notes & Recommendations

* The custom QC filters are designed for **WES (~100–500×)** data
* For low-VAF detection (e.g., AML subclones), consider:

  * lowering AD threshold
  * adding VAF-based filtering
* Ensure normal samples have sufficient coverage to avoid false positives

---

## Future Improvements

* Skip completed cases automatically
* Add VAF-based filtering option
* Support panel of normals (PoN)
* Multi-caller ensemble integration

---