#!/bin/bash

# A script to add variant annotations using vep 

export file=$1
export output1=$2


/opt/compsci/perl/5.24.0/bin/perl  /data/yadavv/ensembl-vep/vep.pl -i  $file -o out_vep.vcf -offline --af --af_1kg --af_esp  --af_exac  --vcf  --force_overwrite --no_stats --fields AF,AFR_AF,EUR_AF,AMR_AF,EAS_AF,AA_AF,EA_AF,ExAC_AC,ExAC_AF -fork 4 --cache --pick_allele  --dir_cache /data/shared/vep-cache
cat out_vep.vcf|grep "#" > head.txt 
echo -e  "##INFO=<ID=ESP6500_EA_AF,Number=A,Type=Float,Description="Field 'ESP6500_EA_AF' "> \n ##INFO=<ID=1000Gp3_AFR_AF,Number=A,Type=Float,Description="Field '1000Gp3_AFR_AF' "> \n ##INFO=<ID=ExAC_AC,Number=A,Type=Integer,Description="Field 'ExAC_AC' from "> \n ##INFO=<ID=1000Gp3_AF,Number=A,Type=Float,Description="Field '1000Gp3_AF' "> \n ##INFO=<ID=ESP6500_AA_AF,Number=A,Type=Float,Description="Field 'ESP6500_AA_AF' "> \n ##INFO=<ID=1000Gp3_AMR_AF,Number=A,Type=Float,Description="Field '1000Gp3_AMR_AF' "> \n ##INFO=<ID=1000Gp3_EUR_AF,Number=A,Type=Float,Description="Field '1000Gp3_EUR_AF' "> \n ##INFO=<ID=ExAC_AF,Number=A,Type=Float,Description="Field 'ExAC_AF' "> \n ##INFO=<ID=1000Gp3_EAS_AF,Number=A,Type=Float,Description="Field '1000Gp3_EAS_AF'">" >> head.txt



cat out_vep.vcf|grep -v "#"|awk -F$'\t' '{ for ( n=1; n<=NF; n++ ) {if($n ~"CSQ=")c1=n;}}NR>0{split($c1,a,"|"); print $1,$2,$3,$4,$5,$6,$7,$8";1000Gp3_AF="a[1]";1000Gp3_AFR_AF="a[2]";1000Gp3_EUR_AF="a[3]";1000Gp3_AMR_AF="a[4]";1000Gp3_EAS_AF="a[5]";ESP6500_AA_AF="a[6]";ESP6500_EA_AF="a[7]";ExAC_AC="a[8]";ExAC_AF="a[9],$9}'|tr ' ' '\t' > content.txt 
cat head.txt content.txt > $output1
rm head.txt content.txt out_vep.vcf