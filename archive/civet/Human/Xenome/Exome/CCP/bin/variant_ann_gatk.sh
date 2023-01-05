#!/bin/bash

# A script to add variant annotations using vep 

export file=$1
export output1=$2

vep -verbose -debug -i  $file -o out_vep.vcf -offline --af --af_1kg --af_esp  \
    --af_exac --vcf  --force_overwrite --no_stats \
    --fields AF,AFR_AF,EUR_AF,AMR_AF,EAS_AF,AA_AF,EA_AF,ExAC_AC,ExAC_AF \
    -fork 4 --cache --pick_allele  --dir_cache /data/shared/vep-cache \
    --cache_version=87

status=$?

if [[ -s out_vep.vcf ]] ; then
    echo "The VEP step ran successfully"
else
    echo "The VEP step crashed with exit status ${status}. Please try a rerun. " >&2
    exit 1
fi ;

<<<<<<< HEAD
cat > $output1 << EOF
=======
cat > output1 << EOF
>>>>>>> f9186f7b6840242cbe1f17521398e82453ee3578
##INFO=<ID=CALLER,Number=1,Type=String,Description="The variant caller Used to call the variant GATK  or Pindel">
##INFO=<ID=ESP6500_EA_AF,Number=A,Type=Float,Description="ESP6500_EA_AF">
##INFO=<ID=1000Gp3_AFR_AF,Number=A,Type=Float,Description="1000Gp3_AFR_AF">
##INFO=<ID=ExAC_AC,Number=A,Type=Float,Description="ExAC_AC">
##INFO=<ID=1000Gp3_AF,Number=A,Type=Float,Description="1000Gp3_AF">
##INFO=<ID=ESP6500_AA_AF,Number=A,Type=Float,Description="ESP6500_AA_AF">
##INFO=<ID=1000Gp3_AMR_AF,Number=A,Type=Float,Description="1000Gp3_AMR_AF">
##INFO=<ID=1000Gp3_EUR_AF,Number=A,Type=Float,Description="1000Gp3_EUR_AF">
##INFO=<ID=ExAC_AF,Number=A,Type=Float,Description="ExAC_AF">
##INFO=<ID=1000Gp3_EAS_AF,Number=A,Type=Float,Description="1000Gp3_EAS_AF">
EOF

<<<<<<< HEAD
grep "#" out_vep.vcf >> $output1

grep -v "#" out_vep.vcf | awk -F$'\t' '{ for ( n=1; n<=NF; n++ ) {if($n ~"CSQ=")c1=n;}}NR>0{split($c1,a,"|"); split(a[2],b,"&"); split(a[3],c,"&"); split(a[4],d,"&"); split(a[5],e,"&"); split(a[6],f,"&"); split(a[7],g,"&"); split(a[8],h,"&"); split(a[9],i,"&"); split(a[10],j,"|");print $1,$2,$3,$4,$5,$6,$7,$8";1000Gp3_AF="b[1]";1000Gp3_AFR_AF="c[1]";1000Gp3_EUR_AF="d[1]";1000Gp3_AMR_AF="e[1]";1000Gp3_EAS_AF="f[1]";ESP6500_AA_AF="g[1]";ESP6500_EA_AF="h[1]";ExAC_AC="i[1]";ExAC_AF="j[10]";CALLER=GATK",$9}'|tr ' ' '\t' >> $output1
=======
grep "#" out_vep.vcf >> output1

grep -v "#" out_vep.vcf | awk -F$'\t' '{ for ( n=1; n<=NF; n++ ) {if($n ~"CSQ=")c1=n;}}NR>0{split($c1,a,"|"); split(a[2],b,"&"); split(a[3],c,"&"); split(a[4],d,"&"); split(a[5],e,"&"); split(a[6],f,"&"); split(a[7],g,"&"); split(a[8],h,"&"); split(a[9],i,"&"); split(a[10],j,"|");print $1,$2,$3,$4,$5,$6,$7,$8";1000Gp3_AF="b[1]";1000Gp3_AFR_AF="c[1]";1000Gp3_EUR_AF="d[1]";1000Gp3_AMR_AF="e[1]";1000Gp3_EAS_AF="f[1]";ESP6500_AA_AF="g[1]";ESP6500_EA_AF="h[1]";ExAC_AC="i[1]";ExAC_AF="j[10]";CALLER=GATK",$9}'|tr ' ' '\t' >> output1
>>>>>>> f9186f7b6840242cbe1f17521398e82453ee3578

rm -f out_vep.vcf
