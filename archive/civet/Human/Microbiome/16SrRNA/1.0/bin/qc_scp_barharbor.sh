WLA='/data/weinstocklab/raw/archive/'
QC_FOLDER=${1?"Folder to transfer QC files from?"};
LOCATION=${2:-'NYGenome/NYGC_Stats_files/Completed_Projects'} # or 'Runs'
QC_ARCH=${WLA}/${QC_FOLDER}_qc/run;
mkdir -p ${QC_ARCH}
scp -r cadillac:/illumina_temp/${LOCATION}/${QC_FOLDER}/{InterOp,RunInfo.xml,runParameters.xml,SampleSheet.csv} ${QC_ARCH};
