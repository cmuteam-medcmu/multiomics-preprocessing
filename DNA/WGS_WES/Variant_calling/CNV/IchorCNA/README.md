# IchorCNA Analysis Pipeline 

## Including tools
-  readCounter
-  IchorCNA

## Setting
### 1. Create/Ensure folder structure (BASE_DIR)
The pipeline expects the output from the Preprocessing Pipeline as follows:

```text
20260204_Batch/ (BASE_DIR)
└── 02_PREPROCESS/
    └── {Sample_ID}/
        └── {Sample_ID}_recal.bam
```

### 2. sample_sheet.csv
    - Use the same CSV file from the Preprocessing step.

| element | format | example |
| :--- | :--- | :--- |
| date | YYYYMMDD | 20260204 |
| library | 3 ABBR. + coverage | WES500X |
| platform | 3 ABBR. | ILL (Illumina) |
| case id | 5 digits number | ST0032 |
| sample type | 3 ABBR. | TSF (Tissue) |
| sample event | disease specific | WES (Whole Exome) |

```csv
ID,R1,R2
20260204_WES500X_ILL_ST0032_TSF_WES,/path/to/R1.fq.gz,/path/to/R2.fq.gz
```

3. shared-IchorCNA.sh
    -  Configure SLURM header for resource allocation.
    -  Set project-specific variables and the path to sample_sheet.csv.
  
```bash
#SBATCH --account=account              
#SBATCH --job-name=IchorCNA           
#SBATCH --partition=compute              
#SBATCH --time=04:00:00                
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G                
#SBATCH --ntasks-per-node=1             
#SBATCH --output=IchorCNA_%j.out         
#SBATCH --error=IchorCNA_%j.err  
```

```bash
readonly BASE_DIR="/project/o260002_CFSTS69/Tissue/WES-500x/20260204_HN00265200_EXO_Report"
readonly SAMPLES_SHEET="${BASE_DIR}/sample_sheet.csv"
readonly GENOME_BUILD="hg38"
```

4. Final output 
'''text
20260204_Batch/ (BASE_DIR)
├── 02_PREPROCESS/
│   └── {Sample_ID}/
│       └── {Sample_ID}_recal.bam          # Input
└── 03_ICHORCNA/
    └── {Sample_ID}/
        ├── {Sample_ID}.bin                # HMMcopy bin-level read counts
        ├── {Sample_ID}.cna.seg            # IchorCNA segmentation results
        ├── {Sample_ID}.seg                # IGV-compatible segment file
        ├── {Sample_ID}.params.txt         # Estimated tumor fraction and ploidy
        ├── {Sample_ID}_genomeWide.pdf     # Genome-wide copy number profile plot
        └── {Sample_ID}_bias.pdf           # GC and mappability bias plots
'''
