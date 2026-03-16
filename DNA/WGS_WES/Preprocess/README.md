# Preprocessing Pipeline 
## Including tools
-  Fastp
-  Fastqc
-  Fq2Bam
-  MosMepth
-  MultiQC

## Setting
1.  sample_sheet.csv
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
