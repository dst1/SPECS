#!/bin/bash

while read old_id new_id
do
	if [ ! -z "$old_id" ]
	then
        	remap_id[${old_id}]=${new_id}
	fi
done < sample_remap.txt

cp raw_data/*.fastq Results/fastq/

cd Results/fastq

for i in "${!remap_id[@]}"
do
	echo "$i --- ${remap_id[$i]}"
	mv ${i}.fastq ${remap_id[$i]}.fastq
done
