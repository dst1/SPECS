#!/bin/bash

s=$1

#first - collapse reads
i=Results/fastq_clip
o=Results/collapsed

echo "collapsing sample ${s}"
mkdir $o
fastx_collapser -Q33 -i ${i}/clipped_${s}.fastq -v -o ${o}/collapsed_${s}.fasta
echo "collapse completed - sample ${s}"

#now align with bowtie
i=Results/collapsed
o=Results/bowtie
ind=Bowtie2_inds

mkdir $o
mkdir $o/Sample_${s}

echo "Sample no. ${s}"
mkdir $o/Sample_${s}/D5
#align reads
echo "Align sample ${s} with bowtie2"
bowtie2 --very-sensitive \
  --norc -f -p 32 \
  --met-file $o/Sample_${s}/D5/met_${s}.txt \
  -x $ind/inds/D5 \
  -U $i/collapsed_${s}.fasta \
  -S $o/Sample_${s}/D5/${s}_D5.SAM 2>&1 | tee $o/Sample_${s}/D5/log_${s}_D5.txt
#filter aligned reads
echo "filtering mapped reads"
samtools view -F 4 $o/Sample_${s}/D5/${s}_D5.SAM | cut -f1,3 > $o/Sample_${s}/D5/mapped_${s}_D5.txt
echo "sample ${s} complete"

