#!/bin/bash

f=$1

echo $f

i=Results/fastq
o=Results/FastQC_pre

mkdir $o
mkdir $o/Sample_$f

fastqc -o "${o}/Sample_${f}" -t 4 --extract "${i}/${f}.fastq"
