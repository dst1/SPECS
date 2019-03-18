# Preprocessing

**All scripts can be ran using:**
`bash RunAll.sh PATH/TO/WORKING/DIR PATH/TO/SCRIPTS/DIR`

The pipeline requires 2 thing:
- raw fastq files in the folder raw_data
- samples file - One sample name per line (same as file name - `SAMP.fastq`) 

<<<<<<< HEAD
## 1. Fastqc
=======
Preprocessing should take approximately 5 minutes per sample (tested on a Dell XPS 9550 laptop running Kubuntu 18.04)

## 1. Renaming
input: fastq files in the folder raw_data
output: fastq files in the folder `Results/fastq` renamed according to remaping file

## 2. Fastqc
>>>>>>> a3ceb112468975a65939a46260b04aa6f65707ba
Fastqc before clipping:
input: fastq files with sample name in the folder `Results/fastq/$samp.fastq`
output: fastqc output in the folder `Results/FastQC_pre/Sample_$samp`

`bash fastqc_pre.sh $samp`

## 2. Clipping
clipping using fastx_clipper and the restriction site, leaving only reads with the restriction site:
input: fastq files with sample name in the folder `Results/fastq`
output: fastq files with the prefix `clipped_$samp` in `Results/fastq_clip`

`bash clipper.sh $samp`

## 3. Fastqc of post clipped reads
Fastqc after clipping:
input: fastq files with the prefix `clipped_$samp` in the folder `Results/fastq_clipped`
output: fastqc output in the folder `Results/FastQC_post/Sample_$samp`

`bash fastqc_post.sh $samp`

## 4. Collapsing reads
collapsing using fastx_collapser
input: fastq files with the prefix `clipped_$samp` in `Results/fastq_clip`
output: fasta files with the prefix `collapsed_$samp` in `Results/collapsed`

## 5. Alignment
Bowtie2 alignment
input: fasta files with the prefix c`ollapsed_$samp` in `Results/collapsed`
output: SAM file all alignments with in the folder `bowtie/Sample_$samp/D5/`

## 6. Filtering
Samtools filter mapped
input: SAM file all alignments in the folder `bowtie/Sample_$samp/D5/`
output: tab delimited file with a line for each seq mapped and its read counts with the prefix `mapped_$samp_D5` in the folder `bowtie/Sample_$samp/D5/`

*(running steps 4-6)*
`bash Bowtie_run.sh $samp`

## 7. Counts table generation
This step is manually executed using the Rmd notebook `Counts_processing.Rmd`
This step also performs counts normalization using DESeq2 and save them to file.

-------------

bowtie indexes for steps 6 and 8 can be found in the dropbox folder:
[Dropbox Bowtie](https://www.dropbox.com/sh/9kbei6qvwkbbtam/AACA3az8v1Ie-wq446Yp8w52a?dl=0)

sample input and output can be found in the dropbox folder:
[Dropbox Example](https://www.dropbox.com/s/0wlx800a9c3nseq/Results.tar.gz?dl=0)