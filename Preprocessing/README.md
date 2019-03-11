# Preprocessing

**All scripts can be ran using:**
`bash RunAll.sh PATH/TO/WORKING/DIR PATH/TO/SCRIPTS/DIR`

The pipeline requires 2 thing:
- raw fastq files in the folder raw_data
- remaping file - first colum contains fastq file names and the second the sample name, space delimited and without headers

## 1. Renaming
input: fastq files in the folder raw_data
output: fastq files in the folder `Results/fastq` renamed according to remaping file

## 2. Fastqc
Fastqc before clipping:
input: fastq files with sample name in the folder `Results/fastq/$samp.fastq`
output: fastqc output in the folder `Results/FastQC_pre/Sample_$samp`

`bash fastqc_pre.sh $samp`

## 3. Clipping
clipping using fastx_clipper and the restriction site, leaving only reads with the restriction site:
input: fastq files with sample name in the folder `Results/fastq`
output: fastq files with the prefix `clipped_$samp` in `Results/fastq_clip`

`bash clipper.sh $samp`

## 4. Fastqc of post clipped reads
Fastqc after clipping:
input: fastq files with the prefix `clipped_$samp` in the folder `Results/fastq_clipped`
output: fastqc output in the folder `Results/FastQC_post/Sample_$samp`

`bash fastqc_post.sh $samp`

## 5. Collapsing reads
collapsing using fastx_collapser
input: fastq files with the prefix `clipped_$samp` in `Results/fastq_clip`
output: fasta files with the prefix `collapsed_$samp` in `Results/collapsed`

## 6. Alignment
Bowtie2 alignment
input: fasta files with the prefix c`ollapsed_$samp` in `Results/collapsed`
output: SAM file all alignments with in the folder `bowtie/Sample_$samp/D5/`

## 7. Filtering
Samtools filter mapped
input: SAM file all alignments in the folder `bowtie/Sample_$samp/D5/`
output: tab delimited file with a line for each seq mapped and its read counts with the prefix `mapped_$samp_D5` in the folder `bowtie/Sample_$samp/D5/`

*(running steps 5-7)*
`bash Bowtie_run.sh $samp`

## 8. Counts table generation
This step is manually executed using the Rmd notebook `Counts_processing.Rmd`
This step also performs counts normalization using DESeq2 and save them to file.

-------------

bowtie indexes for steps 6 and 8 can be found in the dropbox folder:
[Dropbox](https://www.dropbox.com/sh/9kbei6qvwkbbtam/AACA3az8v1Ie-wq446Yp8w52a?dl=0)
