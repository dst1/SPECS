#!/bin/bash

f=$1

i=Results/fastq_clip
o=Results/FastQC_post

mkdir $o
mkdir $o/Sample_$f

fastqc -o "${o}/Sample_${f}" -t 4 --extract "${i}/clipped_${f}.fastq"
