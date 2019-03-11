#!/bin/bash

s=$1
f=$(sed -n ${s}p sample_remap.txt | cut -f2)

i=Results/fastq_clip
o=Results/FastQC_post

mkdir $o
mkdir $o/Sample_$f

fastqc -o "${o}/Sample_${f}" -t 4 --extract "${i}/clipped_${f}.fastq"
