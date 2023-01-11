#!/bin/bash


export file=$1
export output=$2

# formatting Ensembl main table

sed '1,2d' $file | awk '{OFS="\t";$2=$3=$4=$5=$6="";$0=$0;$1=$1}1' | grep -w 'lincRNA\|miRNA\|misc_RNA\|Mt_rRNA\|Mt_tRNA\|snRNA\|snoRNA\|protein_coding\|rRNA\|processed_transcript\|pseudogene\|sense_overlapping' | sort -k 1 > $output
sed -i "1iBioTypes\tExpressionCounts" $output

