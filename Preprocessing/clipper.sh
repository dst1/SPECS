#fastx_clipper

#clipping to keep only the promotor sequences
#removing adapters

s=$1
f=$(sed -n ${s}p sample_remap.txt | cut -f2)

i=Results/fastq
o=Results/fastq_clip

mkdir $o

echo "Processing sample $f..."
# take action on each file. $f store current file name
# Asc1 - GGCGCGC
srun fastx_clipper -a GGCGCGCC -Q33 -c -l 30 -i ${i}/${f}.fastq -o ${o}/clipped_${f}.fastq

