process GATK_FILTERMUTECTCALLS {
  tag "$sampleID"

  cpus 1
  memory 15.GB
  time '05:00:00'

  container '/pdx/pdx_resource_service/elion/containers/gatk-4.0.5.1_htslib_tabix.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", mode:'copy'

  input:
  tuple val(sampleID), file(vcf), file(tbi)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf
  file("*.vcf.idx")

  script:
  
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

    """
    java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar /gatk/gatk.jar \
    FilterMutectCalls \
    --variant ${sampleID}_intermed.vcf.gz \
    --output ${sampleID}_final_mutect_snp_indel_filtered.vcf \
    --min-base-quality-score 20 \
    --dont-use-soft-clipped-bases true \
    --unique-alt-read-count 5 \
    --stand-call-conf 30 \
    -L ${params.targets_gatk}
    """
}
