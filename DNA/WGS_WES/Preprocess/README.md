# Preprocessing Pipeline 
## Including tools
-  Fastp
-  Fastqc
-  Fq2Bam
-  MosMepth
-  MultiQC

## Setting
1. sample_sheet.csv
    -  ID must follow AWS structure in each project

  |element|format|example|
  |----|----|----|
  |date|YYYYMMDD|20260101|
  |library|3 ABBR. with/without coverage|WGS20X|
  |platform|3 ABBR.|ILU|
  |case id|5 digits number|OS00001|
  |sample type|3 ABBR.|CFD|
  |sample event|based on each disease|Fx12M|

```csv
ID,R1,R2
20260101_WGS20X_ILU_OS00001_CFD_Fx12M,/path/to/Illumina_R1.fq.gz,/path/to/Illumina_R2.fq.gz
```

2. run_preprocess_pipeline.csv
    -  Change header for slurm management
    -  Set variables and the paths of sample_sheet.csv & qc_script.py
  
```bash
#SBATCH --account=o250022               
#SBATCH --job-name=fq2bam  
#SBATCH --partition=gpu
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=40
#SBATCH --gpus=1
#SBATCH --mem=200G
#SBATCH --output=preprocess_%j.out          
#SBATCH --error=preprocess_%j.err
```

```bash
readonly THREADS=40
readonly MEMORY=200 # Gb unit
readonly OUT_DIR="/project/o250022_cfOSTEO/script/preprocess"
readonly SAMPLES_SHEET="./sample_sheet.csv"
readonly QC_SCRIPT="./qc_script.py"
readonly REPORT_DIR="${OUT_DIR}/reports"
```
   
3. QC template
    - Copy data from QC_SUMMARY.csv
    - Paste to Google sheet template [QC Summary](https://docs.google.com/spreadsheets/d/1JUirnX6yo0NkUpjTHGlRnbIpUlB-UyOtZFrjM2WKcJs/edit?usp=sharing).
