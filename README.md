# SPECS
Code repository for the article: A High-Throughput Screening and Computation Platform for Identifying Synthetic Promoters with Enhanced Cell- State Specificity (SPECS)

# Structure
- **Preprocessing** - contains scripts for the preprocessing of SPECS experiments, starting from FASTQ files and ending with normalized counts

- **Fluorescence prediction** - contains a Rmarkdown notebook that used to process FACS readings together with normalized counts to produce a machine learning predictive model of fluorescence.

# Requirements

For pre-processing:
```
bowtie2
fastx_toolkit
fastqc
SAMtools
```

For the predictive model:

```
R 3.4 or greater

R Packages:
-----------

DESeq2
data.table
reshape2
tidyverse
flowCore
caret
```
