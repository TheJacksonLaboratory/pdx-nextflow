 
#!/bin/bash

# A script to Add Caller column

export file=$1
export output1=$2




cat $file|grep "#" > head.txt 
cat $file|grep "GT:AD:DP:GQ:PL" > tail.txt 

cat head.txt tail.txt > $output1

rm head.txt tail.txt
