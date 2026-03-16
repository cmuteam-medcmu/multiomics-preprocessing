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
  -  {date in yyyymmdd}_{platform with coverage}_{case id with 5 digits}_{sample type}_{sample event}
```csv
ID,R1,R2
20260101_WGS20X_OS00007_CF_Fx12M,/path/to/Illumina_R1.fq.gz,/path/to/Illumina_R2.fq.gz
```
