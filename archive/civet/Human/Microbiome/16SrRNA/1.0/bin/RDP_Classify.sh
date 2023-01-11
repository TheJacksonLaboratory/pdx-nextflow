#!/usr/bin/env bash
# QSUB Vars
#PBS -N RDP_Classify
#PBS -M benjamin.leopold@jax.org
#PBS -m abe
#PBS -j oe
#PBS -q batch
#PBS -l nodes=1:ppn=1
#PBS -l mem=8GB
#PBS -l walltime=8:00:00

module add fastx

# bin Variables
CLASSIFIER_JAR_PATH="/opt/compsci/Microbiome/16srRNA/1.0/RDPTools/classifier.jar"
FASTX_PATH="fastq_to_fasta"

#Script Variables
BASE_DIR=${BASE_DIR:-"/data/pbais/Microbiome/Blake/TrialRunBatch"}
#CLASS_DIR="${BASE_DIR}/${2:-"classified"}"
FASTA_DIR="${BASE_DIR}/${3:-"fasta"}"

CONVERT_FASTQ=${CONVERT_FASTQ:-"yes, please"}
#CONVERT_FASTQ=""

#FASTQ_FILE=${1:?"pass fastq file as first arg!"}
FASTQ_FILE=${FASTQ_FILE:?"set env var!"}
FASTA_FILE=$( basename ${FASTQ_FILE/%q/a} )

#LOGFILE=RDPc_${FASTA_FILE%.fasta}.log
CLASSIFIED_OUT_FILE=${CLASS_DIR}/${FASTA_FILE%.fasta}_classified.tsv
HIER_OUTFILE=${CLASS_DIR}/${FASTA_FILE%.fasta}_taxa.tsv

# RDP vars
CUTOFF=${CUTOFF:-'0.5'}
FORMAT=${FORMAT:-'filterbyconf'}

#Code
cd $BASE_DIR
mkdir -p $CLASS_DIR $FASTA_DIR

# format convert:
if [ -n "$CONVERT_FASTQ" ] ; then
    echo "Converting fastq to fasta format: $FASTQ_FILE"
    $FASTX_PATH -Q33 -i $FASTQ_FILE -o $FASTA_DIR/$FASTA_FILE
fi

#Note: -c and -f default values are '0.8' and 'allrank'
echo "Classifying $FASTA_FILE using RDP classifier using cutoff of $CUTOFF"
java -Xmx1g -jar $CLASSIFIER_JAR_PATH classify \
        -c $CUTOFF \
        -o $CLASSIFIED_OUT_FILE \
        -f $FORMAT \
        -h $HIER_OUTFILE \
        $FASTA_DIR/$FASTA_FILE

#        2>&1 >> $LOGFILE

# convert TSV to CSV
echo "Converting RDP results from Tab to Comma delimiters"
tr "\t" "," < $CLASSIFIED_OUT_FILE > ${CLASSIFIED_OUT_FILE%.tsv}.csv
tr "\t" "," < $HIER_OUTFILE > ${HIER_OUTFILE%.tsv}.csv
rm $CLASSIFIED_OUT_FILE $HIER_OUTFILE

echo # empty spacer line at end of run
