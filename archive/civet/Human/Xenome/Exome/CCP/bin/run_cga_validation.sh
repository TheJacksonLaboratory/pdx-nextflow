#! /bin/bash

#
# Run the three cga validation test sets for single_sample_exome
#
# NOTE!  Must be invoked with the version number of cga that is to be tested, e.g.,
#
#     $ run_cga_validation.sh 1.2.2
#
#

if [ "X${1}" = "X" ]; then
    echo "USAGE: ${0} <cga-version-number>"
    exit 1
fi

if [ ! -d /opt/compsci/cga/${1} ]; then
    echo "ERROR: ${1} does not appear to be a valid cga version number."
    echo "       could not find /opt/compsci/cga/${1}/"
    exit 2
fi

module unload cga civet
module load cga/${1}
module list
echo "Testing cga and civet with path:"
echo $PATH

FQD=/sequenced-data/clinical/cga/validation/single_sample_exome/fastq

DIR=testing_results/${1}/uncompressed
echo "Testing uncompressed: results in ${DIR}"
mkdir -p ${DIR}
pushd ${DIR}
pwd
single_sample_exome ${FQD}/NA12892NA18507_S1_L001_R{1,2}_001.fastq
single_sample_exome ${FQD}/NA12892NA18507_S2_L001_R{1,2}_001.fastq
single_sample_exome ${FQD}/NA12892NA18507_S3_L001_R{1,2}_001.fastq
popd

# Now test fastq.gz processing
DIR=testing_results/${1}/compressed
echo ""
echo "Testing compressed: results in ${DIR}"
mkdir -p ${DIR}
pushd ${DIR}
pwd 
single_sample_exome ${FQD}/NA12892NA18507_S1_L001_R{1,2}_001.fastq.gz
single_sample_exome ${FQD}/NA12892NA18507_S2_L001_R{1,2}_001.fastq.gz
single_sample_exome ${FQD}/NA12892NA18507_S3_L001_R{1,2}_001.fastq.gz
popd
