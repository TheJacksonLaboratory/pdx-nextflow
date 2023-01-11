#!/bin/bash


export file=$1
export output=$2

fbname=$(basename "$file" _quantification.txt)

sed '1,2d' $file | awk '{OFS="\t";$2=$3=$4=$5=$6="";$0=$0;$1=$1}1' > $output
sed -i "1imiRNA\t$fbname" $output

##sed '1,2d' {out_1} | awk '{OFS="\t";$2=$3=$4=$5=$6="";$0=$0;$1=$1}1' > {out_2}
