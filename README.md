# SPECS
Code repository for the article: A High-Throughput Screening and Computation Platform for Identifying Synthetic Promoters with Enhanced Cell- State Specificity (SPECS) https://www.nature.com/articles/s41467-019-10912-8

citation (_bibtex_):
```{bibtex}
@article{Wu2019,
  doi = {10.1038/s41467-019-10912-8},
  url = {https://doi.org/10.1038/s41467-019-10912-8},
  year = {2019},
  month = jun,
  publisher = {Springer Science and Business Media {LLC}},
  volume = {10},
  number = {1},
  author = {Ming-Ru Wu and Lior Nissim and Doron Stupp and Erez Pery and Adina Binder-Nissim and Karen Weisinger and Casper Enghuus and Sebastian R. Palacios and Melissa Humphrey and Zhizhuo Zhang and Eva Maria Novoa and Manolis Kellis and Ron Weiss and Samuel D. Rabkin and Yuval Tabach and Timothy K. Lu},
  title = {A high-throughput screening and computation platform for identifying synthetic promoters with enhanced cell-state specificity ({SPECS})},
  journal = {Nature Communications}
}
```


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
