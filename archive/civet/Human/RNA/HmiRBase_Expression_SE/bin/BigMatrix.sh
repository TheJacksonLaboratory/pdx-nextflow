#!bin/bash


## Final formating to create entire expression into one matrix

export BigMatrix=$1
export first=$2
export second=$3
export output=$4

find . -type f -name '*_miR_Counts.txt' -exec paste {} + > $BigMatrix

colmn=$(awk '{print NF}' $BigMatrix | tail -n 1)

if [[ "$colmn" -gt 2 ]]; then

	awk '{print $1}' $BigMatrix > $first

	awk '{sept=""; for (i=2;i<=NF;i+=2) {printf "%s%s", sept, $i; sept="\t"}; printf "\n" }' $BigMatrix > $second

	paste $first $second > $output

	rm $BigMatrix $first $second
else
	rm $BigMatrix
fi

