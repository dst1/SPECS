#!/bin/bash

work_dir=$1
scripts_dir=$2

cd ${work_dir}

mkdir -p Results
mkdir -p Results/fastq

for samp in $(cat samples.txt)
do
    echo $samp
    #FastQC Pre:
    mkdir -p Results/FastQC_pre
    bash ${scripts_dir}/fastqc_pre.sh $samp

    #Clipping:
    mkdir -p Results/fastq_clip
    bash ${scripts_dir}/clipper.sh $samp

    #FastQC Post:
    mkdir -p Results/FastQC_post
    bash ${scripts_dir}/fastqc_post.sh $samp

    #Collapse, align, filter:
    mkdir -p Results/collapsed
    mkdir -p Results/bowtie
    
    bash ${scripts_dir}/Bowtie_run.sh $samp
done
