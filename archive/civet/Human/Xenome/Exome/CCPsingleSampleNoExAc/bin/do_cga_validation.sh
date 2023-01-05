#! /bin/bash

# This script runs the comparison phase of the CGA validation 
# process specified in SOP.CBC.005 Testing CGA Pipeline.docx.
# It should produce no output.

# Invoke the script with a single parameter, the prefix to the
# output directories, through the yymmdd.  For example:
#     do_cga_validation.sh 130822
# or
#     do_cga_validation.sh /home/asimons/projects/cga/130822

# This script should produce no output

# Check for args:
if [ "X${2}" = "X" ]; then
    echo "USAGE: ${0} <cga-version> <yymmdd-date>"
    exit 1
fi

if [ ! -d /opt/compsci/cga/${1} ]; then
    echo "ERROR: ${1} does not appear to be a valid cga version number."
    echo "       could not find /opt/compsci/cga/${1}/"
    exit 2
fi

BM=/sequenced-data/clinical/cga/validation/single_sample_exome/benchmark_results
FILE1=NA12892NA18507_S1_L001
FILE2=NA12892NA18507_S2_L001
FILE3=NA12892NA18507_S3_L001
END1=variants_microIndels.DPfiltered.clinicalTargets.vcf
END2=variants_microIndels.DPfiltered.clinicalTargets.AFge0.05.tab
END3=variants_microIndels.DPfiltered.Annotated.tab

# This file is in CNV_output
END4=CNATable.50rd.20bases.20bins.txt

VALUES=stat_qual_summary_values.txt

echo "This test should only output two lines starting with \"Checking\""
echo "Any other output after the \"START TEST\" line indicates a failure."
echo "START TEST"

# The directory for this version's test results.
BASE=testing_results/${1}

for D in compressed uncompressed; do
    B=${BASE}/${D}
    echo "Checking ${B}"
    for F in $FILE1 $FILE2 $FILE3; do
        # The directory for this sample's run:
        SAMP_OUT=${B}/${2}_${F}
        for I in 1 2 3 4; do

            # The following two lines are serious bash voodoo,
            # taken from StackOverflow.
            TE=END${I}
            E=${!TE}  # The ! forces a kind of indirection

            # This mess of having FN and FN_U is because one of our
            # test files is down a directory.  FN has the extra
            # directory; FN_U doesn't.
            FN=${F}_${E}
            FN_U=${F}_${E}
            if [ ${I} -eq 4 ]; then
                FN=CNV_output/${F}.${E}
                FN_U=${F}.${E}
            fi
            LF=${SAMP_OUT}/${FN}
            LFF=TEST_${FN_U}_filtered
            BF=${BM}/${FN_U}
            BFF=BENCHMARK_${FN_U}_filtered
            if [ -e ${LF} ] ; then
                filter_validation_vcf.py ${LF} > ${LFF}
                filter_validation_vcf.py ${BF} > ${BFF}
                diff ${BFF} ${LFF}
                rm ${BFF} ${LFF}
            else
                echo "FAILURE: ${LF} does not appear to exist."
            fi
        done

        # Benchmark values file
        BV=${BM}/${F}_${VALUES}
        BVF=BENCHMARK_${F}_${VALUES}_filtered
        # Test values file
        LV=${SAMP_OUT}/${VALUES}
        LVF=TEST_${F}_${VALUES}_filtered

        if [ -e ${LV} ] ; then
            # Trim off columns we don't want to compare
            cut -f 3-23 ${BV} > ${BVF}
            cut -f 3-23 ${LV} > ${LVF}
            diff ${BVF} ${LVF}
            rm ${BVF} ${LVF}
        else
            echo "FAILURE: ${LF} does not appear to exist."
        fi
    done
done

RUN_CHECK=${BASE}/check_run.txt
RUN_VALUES=${BASE}/check_run_values.txt
/opt/compsci/cga/${1}/bin/check_cga_run.py -o ${RUN_CHECK} ${BASE} >/dev/null 2>&1
if ! ls ${RUN_CHECK} > /dev/null 2>&1
then
    echo "Could not find ${RUN_CHECK}"
else
    grep "ERROR\|FAIL" ${RUN_CHECK}
fi
if ! ls ${RUN_VALUES} > /dev/null 2>&1
then
    echo "Could not find ${RUN_VALUES}"
else
    grep "ERROR\|FAIL" ${RUN_VALUES}
fi

