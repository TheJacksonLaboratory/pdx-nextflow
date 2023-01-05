#!/bin/bash

# A script to format sort VCF file  according to reference to avoid errors in downstream 

export file=$1
export output=$2
module load tabix
export chrlist=/data/shared/cga_reference_data/hg38_201601/hg38chr_list.txt
bgzip $file
tabix -pvcf $file.gz
cat $chrlist|xargs tabix -h  $file.gz > $output 


