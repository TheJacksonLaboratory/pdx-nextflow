#!/usr/bin/env bash

USAGE="$0 sample_sheet_file_path fastq_dir"

SAMPLE_SHEET=${1?"
$USAGE"}
FASTQ_DIR=${2?"
$USAGE"}

if [[ -r $SAMPLE_SHEET ]]; then
    echo "Reading sample sheet: $SAMPLE_SHEET"
    SAMPLE_RECORDS=( $(grep -A 150 '^Sample_ID' $SAMPLE_SHEET | grep -v 'Sample_ID') );
    FASTQ_FILES=$( ls ${FASTQ_DIR}/*[._]R1[._]*.fastq.gz | grep -v [Uu]ndetermined );
    NEXTSEQ_DIR_RE="[0-9][0-9][0-9][0-9][0-9][0-9]_NS[0-9][0-9][0-9][0-9][0-9][0-9]_";
    if [[ "$SAMPLE_SHEET" =~ "${NEXTSEQ_DIR_RE}" ]]; then
        FASTQ_FILES=$( ls ${FASTQ_DIR}/*_[SL][0-9]*_R1[._]*.fastq* | grep -v [Uu]ndetermined );
    fi;

    NUM_SAMPLES="${#SAMPLE_RECORDS[*]}"
    echo "SampleSheet lines: $NUM_SAMPLES"

    NUM_FILES="${#FASTQ_FILES[*]}"
    echo "Number of FASTQ files: $NUM_FILES"

    if [[ $NUM_SAMPLES == $NUM_FILES ]]; then
        echo "Number of sample files matches sample sheet."
        IFS=$'\n' SAMPLES=($(sort <<<"${SAMPLE_RECORDS[*]}"));
        IFS=$'\n' FASTQS=($(sort <<<"${FASTQ_FILES[*]}"));
        for (( i = 0; i < ${#SAMPLES}; i++ )); do
            if [[ ! "${FASTQS[$i]}" == "${SAMPLES[$i]}"* ]]; then
                echo 'ERROR: Mismatch between SampleSheet and Fastq file!!';
                echo -e "\tSample: '${SAMPLES[$i]}' and file: '${FASTQS[$i]}'\n" >&2;
                echo -e "\t(Please check this for errors.)\n" >&2;
                exit 11;
            else
                echo "SampleSheet matches number and names of Fastq files! Hooray!" >&1;
                exit 0;
            fi;
        done;
    else
        echo "Uh-Oh!  Number of sample files does NOT match the number in the SampleSheet!!"
    fi
else
    echo "ERROR!!! SampleSheet is not found... '"$SAMPLE_SHEET"'" >&1;
    exit 0;
fi

