process GET_MODEL_GENDER {

    tag "$sampleID"

    cpus 1
    memory 6.GB
    time 8.h
    errorStrategy 'finish'

    container '/pdx/pdx_resource_service/elion/containers/apt2.11.3_python2.7.11.sif'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'apt' }", pattern: "*.{txt,log}", mode:'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'apt' }", pattern: "listfile1", mode:'copy'

    input:
        tuple val(sampleID), path(cel_file)

    output:
        tuple val(sampleID), file("listfile1"), emit: cel_list
        tuple val(sampleID), file("gender.txt"), emit: gender
        file("birdseed.calls1.txt")
        file("birdseed.confidences1.txt")
        file("birdseed.report1.txt")
        file("gender_birdseed.txt")
        file("gender_elims.txt")
        file("*log") //this might not be needed. Output was not originally being saved. 

    script:
    """
    sname=\$(echo ${sampleID} | cut -d"_" -f1)

    python ${projectDir}/bin/cnv/get_model_gender.py \${sname} > gender_elims.txt
 
    echo "cel_files" > listfile

    echo $cel_file >> listfile

    awk -v a=${params.hapmap_dat} '{ if (NR>1) print a"/"\$2}' ${params.hapmap_fm} >> listfile

    awk -v a=${params.hapmap_dat} '{ if (NR>1) print a"/"\$2}' ${params.hapmap_m} >> listfile

    /apt_2.11.3_linux_64_bit_x86_binaries/bin/apt-probeset-genotype -c ${params.snp6chip} -a birdseed --read-models-birdseed ${params.snp6chip_birdseed_mod} --special-snps ${params.snp6chip_specsnps} --out-dir .  --cel-files listfile

    cat birdseed.report.txt | grep -v "#" | awk 'NR==2' | cut -f2 > gender_birdseed.txt

    gender=\$(cat gender_elims.txt)

    if [ "X\$gender" = "Xunknown" -o "X\$gender" = "Xunspecified" ]; then gender=\$(cat gender_birdseed.txt); cp gender_birdseed.txt gender.txt; else cp gender_elims.txt gender.txt; fi

    echo "cel_files" > listfile1

    echo $cel_file >> listfile1

    if [ "X\$gender" = "Xfemale" -o "X\$gender" = "Xunknown" ]; then awk -v a=${params.hapmap_dat} '{if (NR>1) print a"/"\$2}' ${params.hapmap_fm} >> listfile1; elif [ "X\$gender" = "Xmale" ]; then awk -v a=${params.hapmap_dat} '{if (NR>1) print a"/"\$2}' ${params.hapmap_m} >> listfile1; fi

    cat birdseed.confidences.txt | grep -v "#" | cut -f1-2 > birdseed.confidences1.txt

    cat birdseed.calls.txt | grep -v "#" | cut -f1-2 > birdseed.calls1.txt

    cat birdseed.report.txt | grep -v "#" | head -2 > birdseed.report1.txt

    """

    stub:
    """
    touch listfile1
    touch gender.txt
    touch test.log
    touch birdseed.calls1.txt
    touch birdseed.confidences1.txt
    touch birdseed.report1.txt
    touch gender_birdseed.txt
    touch gender_elims.txt
    """
}

//    rm -rf listfile birdseed.report.txt birdseed.confidences.txt birdseed.calls.txt 
//    This was originally included but isn't needed in a nextflow context.