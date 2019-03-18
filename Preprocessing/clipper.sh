#fastx_clipper

#clipping to keep only the promotor sequences
#removing adapters

f=$1

i=Results/fastq
o=Results/fastq_clip

mkdir $o

echo "Processing sample $f..."
# take action on each file. $f store current file name
# Asc1 - GGCGCGC
fastx_clipper -a GGCGCGCC -Q33 -c -l 30 -i ${i}/${f}.fastq -o ${o}/clipped_${f}.fastq

