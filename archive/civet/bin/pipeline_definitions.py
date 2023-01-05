# pipes = {pipe_name:{fullname:, dir:, pipe_file:, def_option_file:,
#                     args:, allow_groups:, trailing_args:,
#                     help_prolog:, help_epilog:}
# The args field describes the arguments. The len(args) specifies the
# number of required arguments, or, if allow_groups is True, the
# multiples of that many, e.g., single end, 1, paired, 2, TN, 4.
#
# The allow_groups should be False, UNLESS
#   1) The arguments are solely fastq files, or there is a trailing_args
#      key, and
#   2) The pipeline has been adapted to use the filelist tag as
#      the second parameter.

# Our different argument patterns

fqse = ['fastq']
fqpe = ['fastq1', 'fastq2']
fqpeg = ['fastq1_e1', 'fastq1_e2', '[fastq2_e1 fastq2_e2] ...']
fqtn = ['tumor_fq1', 'tumor_fq2', 'norm_fq1', 'norm_fq2']
fqclosest = ['fastq1_e1', 'fastq2_e2', '[fastq2_e1 fastq2_e2] ...',
             'known-SNPs', 'known-Indels']
chipc = ['sample_fq1', 'sample_fq2', 'control_fq1', 'control_fq2']
fqgen = ['fastq', 'genome-list']
fqpegen = ['fastq1', 'fastq2', 'genome-list']
cel_single = ['cel']
cel_paired = ['cel_tumor', 'cel_normal']
microbiome_dir = ['in_dir', 'out_dir']

# for these argument patterns, the generation and sex are passed as a period separated 
# string e.g. "G23.F" for genration 23 female
fqpe_gensex = ['fastq1', 'fastq2', 'GenSex']
fqse_gensex = ['fastq1', 'GenSex']

cel_single = ['cel']
cel_paired = ['cel_tumor', 'cel_normal']
microbiome_dir = ['in_dir','out_dir']

# our dummy filenames for testing...
testing_files = {
    'fqse':  ['single_end_R1.fastq.gz'],
    'fqpe':  ['paired_end_R1.fastq.gz',
              'paired_end_R2.fastq.gz'],
    'fqpeg': ['paired_end_f1_R1.fastq.gz',
              'paired_end_f1_R2.fastq.gz',
              'paired_end_f2_R1.fastq.gz',
              'paired_end_f2_R2.fastq.gz'],
    'fqtn':  ['paired_end_tumor_R1.fastq.gz',
              'paired_end_tumor_R2.fastq.gz',
              'paired_end_normal_R1.fastq.gz',
              'paired_end_normal_R2.fastq.gz'],
    'fqclosest': ['paired_end_f1_R1.fastq.gz',
                  'paired_end_f1_R2.fastq.gz',
                  'paired_end_f2_R1.fastq.gz',
                  'paired_end_f2_R2.fastq.gz',
                  'snps.vcf',
                  'indels.vcf'],
    'chipc': ['sample_fq1.fastq.gz',
              'sample_fq2.fastq.gz',
              'control_fq1.fastq.gz',
              'control_fq2.fastq.gz'],
    'fqpegen': ['paired_end_f1_R1.fastq.gz',
              'paired_end_R2.fastq.gz','BALB_cJ'],
    'fqgen': ['single_end_R1.fastq.gz',
                'BALB_cJ'],
    'cel_single': ['cel_file.CEL'],
    'cel_paired': ['cel_tumor_file.CEL',
                   'cel_normal_file.CEL'],
    'microbiome_dir': ['microbiome_in_dir',
                       'microbiome_out_dir'],
    'fqse_gensex': ['single_end_R1.fastq.gz','G10.M'],

    'fqpe_gensex': ['paired_end_f1_R1.fastq.gz',
                    'paired_end_R2.fastq.gz',
                    'G10.M'],
}


def get_cga_exome_dir():
    cga_ver = None
    for line in open('/opt/compsci/modulefiles/cga/.version'):
        if 'ModulesVersion' in line:
            parts = line.split('"')
            cga_ver = parts[1]
            break
    if not cga_ver:
        raise ValueError('Could not find current CGA version')
    return '/opt/compsci/cga/{}/exome'.format(cga_ver)


pipes = {
    'hwes': {
        'fullname': 'Whole Exome Paired End data from patient only',
        'dir': 'Human/Exome/WholeExomeSingleSample',
        'pipeline': 'WholeExomeSingleSample.xml',
        'option': 'config_file_Whole_Exome_Single_Sample',
        'args': fqpeg,
        'num_in_group': 2,
        'allow_groups': True,
        'testing_files': 'fqpeg',
        'help_prolog': """Whole Exome Paired End data from patient only(hg38 version)""",
    },
    'hwes_en': {
        'fullname': 'Whole Exome Paired End data Engrafted tumor',
        'dir': 'Human/Xenome/Exome/WholeExome',
        'pipeline': 'XenomeWholeExomeSingleSample.xml',
        'option': 'config_file_Whole_Exome_Single_Sample',
        'num_in_group': 2,
        'args': fqpeg,
        'allow_groups': True,
        'testing_files': 'fqpeg',
        'help_prolog': """Whole Exome Paired End data Engrafted tumor(hg38 version)""",
    },
    'hwes_tn': {
        'fullname': 'Whole Exome Paired End data Tumor-Normal',
        'dir': 'Human/Exome/WholeExomeTumorNormal',
        'pipeline': 'WholeExomeTumorNormal.xml',
        'option': 'config_file_tumor_normal_Whole_exome',
        'args': fqtn,
        'allow_groups': False,
        'testing_files': 'fqtn',
        'help_prolog': """Whole Exome Paired End data Tumor-Normal(hg38 version)""",
    },
    'hwes_en_tn': {
        'fullname': 'Whole Exome Paired End data Engrafted Tumor Tumor-Normal',
        'dir': 'Human/Xenome/Exome/WholeExomeTumorNormal',
        'pipeline': 'WholeExomeTumorNormal.xml',
        'option': 'config_file_tumor_normal_Whole_exome',
        'args': fqtn,
        'allow_groups': False,
        'testing_files': 'fqtn',
        'help_prolog': """Whole Exome Paired End data Engrafted Tumor Tumor-Normal(hg38 version)""",
    },
    'hchip': {
        'fullname': 'ChIP-seq Paired End Experimental-Control',
        'dir': 'Human/ChIPseq/ChIPseqPairedEnd',
        'pipeline': 'experimentalControlChIP.xml',
        'option': None,
        'args': chipc,
        'allow_groups': False,
        'testing_files': 'chipc',
        'help_prolog': """ChIP-seq Paired End Experimental-Control(hg38 version)""",
    },
    'hctp': {
        'fullname': 'CTP for Paired End data from patient only',
        'dir': None,  # We will use a function for this
        'pipeline': 'single_sample_exome.xml',
        'option': None,
        'args': fqpeg,
        'num_in_group': 2,
        'allow_groups': True,
        'dir_function': get_cga_exome_dir,
        'testing_files': 'fqpeg',
        'use_old_style_list_param': True
    },
    'hctp_en': {
        'fullname': 'CTP for Paired end data Engrafted Tumor',
        'dir': 'Human/Xenome/Exome/CCP',
        'pipeline': 'XenomeSingleSampleCCP.xml',
        'option': 'config_file_Xenome_CCP_Single_Sample',
        'args': fqpeg,
        'num_in_group': 2,
        'allow_groups': True,
        'testing_files': 'fqpeg',
        'help_prolog': """CTP for Paired end data Engrafted Tumor(hg38 version)""",
    },
    'hctp_en_mutect': {
        'fullname': 'CTP for Paired end data Engrafted Tumor',
        'dir': 'Human/Xenome/Exome/CCPsingleSampleNoExAc',
        'pipeline': 'XenomeSingleSampleCCP_mutect2_NoExAC.xml',
        'option': 'config_file_Xenome_CCP_Single_Sample',
        'args': fqpeg,
        'num_in_group': 2,
        'allow_groups': True,
        'testing_files': 'fqpeg',
        'help_prolog': """CTP for Paired end data Engrafted Tumor(hg38 version) with Mutect2 variant calling without germline events.""",
    },
    'hctp_tn': {
        'fullname': 'CTP for Paired end data Tumor-Normal',
        'dir': 'Human/Exome/CCPTumorNormal',
        'pipeline': 'CCPTumorNormal.xml',
        'option': 'config_file_tumor_normal_Whole_exome',
        'args': fqtn,
        'allow_groups': False,
        'testing_files': 'fqtn',
    },
    'hctp_en_tn': {
        'fullname': 'CTP for Paired end data Engrafted Tumor Tumor-Normal',
        'dir': 'Human/Xenome/Exome/CCPTumorNormal',
        'pipeline': 'CCPTumorNormal.xml',
        'option': 'config_file_CCPTumorNormal',
        'args': fqtn,
        'allow_groups': False,
        'testing_files': 'fqtn',
        'help_prolog': """CTP for Paired end data Engrafted Tumor Tumor-Normal(hg38 version)""",
    },
    'hwgs': {
        'fullname': 'Whole Genome Sequencing from Patient Only',
        'dir': 'Human/WholeGenome/WholeGenomeSingleSample',
        'pipeline': 'WholeGenomeSingleSample.xml',
        'option': 'config_file_whole_genome_single_sample',
        'args': fqpeg,
        'num_in_group': 2,
        'allow_groups': True,
        'testing_files': 'fqpeg',
    },
    'hwgs_en': {
        'fullname': 'Whole Genome Sequencing Engrafted tumor',
        'dir': 'Human/Xenome/WholeGenome/WholeGenomeSingleSample',
        'pipeline': 'XenomeWholeGenomeSingleSample.xml',
        'option': 'config_file_whole_genome_single_sample',
        'args': fqpeg,
        'num_in_group': 2,
        'allow_groups': True,
        'testing_files': 'fqpeg',
        'help_prolog': """Whole Genome Sequencing Engrafted tumor(hg38 version)""",
    },
    'hwgs_tn': {
        'fullname': 'Whole Genome Sequencing Tumor-Normal',
        'dir': 'Human/WholeGenome/WholeGenomeTumorNormal',
        'pipeline': 'WholeGenomeTumorNormal.xml',
        'option': 'config_file_tumor_normal_Whole_genome',
        'args': fqtn,
        'allow_groups': False,
        'testing_files': 'fqtn',
    },
    'hrpe': {
        'fullname': 'RNA Paired End Patient only',
        'dir': 'Human/RNA/RNA_Expression_Estimation_Single_Sample_PE',
        'pipeline': 'RNASeqSingleSamplePE.xml',
        'option':
            'config_file_RSEM_RNA_SEQ_Single_Sample_Expression_Estimation',
        'num_in_group': 2,
        'args': fqpeg,
        'allow_groups': True,
        'testing_files': 'fqpeg',
        'help_prolog': """RNA Paired End Patient only(hg38 version)""",
    },
    'hrpe_en': {
        'fullname': 'RNA Paired End Engrafted tumor',
        'dir': 'Human/Xenome/RNA/RNA_Expression_Estimation_Single_Sample_PE',
        'pipeline': 'XenomeRnaSeqSingleSamplePE.xml',
        'option':
            'config_file_RSEM_RNA_SEQ_Single_Sample_Expression_Estimation',
        'num_in_group': 2,
        'args': fqpeg,
        'allow_groups': True,
        'testing_files': 'fqpeg',
        'help_prolog': """RNA Paired End Engrafted tumor(hg38 version)""",
    },
    'hrpe_dfci_en': {
        'fullname': 'RNA Paired End Engrafted tumor for DFCI capture beds',
        'dir': 'Human/Xenome/RNA/RNA_Expression_Estimation_Single_Sample_PE',
        'pipeline': 'XenomeRnaSeqSingleSamplePE.xml',
        'option': 'config_file_DFCI_RSEM_RNA_SEQ_Single_Sample_Expression_Estimation',
        'num_in_group': 2,
        'args': fqpeg,
        'allow_groups': True,
        'testing_files': 'fqpeg',
        'help_prolog': """RNA Paired End Engrafted tumor for DFCI capture beds(hg38 version)""",
    },
    'hrpe_st_en': {
        'fullname': 'RNA Paired End Engrafted tumor with strandness',
        'dir': 'Human/Xenome/RNA/RNA_Expression_Estimation_Single_Sample_PE',
        'pipeline': 'XenomeRnaSeqSingleSamplePE.xml',
        'option':
            'config_file_RSEM_RNA_SEQ_Strandness_Single_Sample_Expression_Estimation',
        'num_in_group': 2,
        'args': fqpeg,
        'allow_groups': True,
        'testing_files': 'fqpeg',
        'help_prolog': """RNA Paired End Engrafted tumor with strandness(hg38 version)""",
    },

    'hrpe_vt_en': {
        'fullname': 'RNA Paired End Engrafted tumor with variant calling using STAR aligner',
        'dir': 'Human/Xenome/RNA/RNA_Expression_Variant_Calling_STAR_Single_Sample_PE',
        'pipeline': 'XenomeRnaSeqSingleSamplePE_Star_VariantCalling.xml',
        'option':
            'config_file_Star_non_stranded',
        'num_in_group': 2,
        'args': fqpeg,
        'allow_groups': True,
        'testing_files': 'fqpeg',
        'help_prolog': """RNA Paired End Engrafted tumor with variant calling using star aligner and GATK(hg38 version)""",
    },

    'hrse': {
        'fullname': 'RNA Single End Patient Only',
        'dir': 'Human/RNA/RNA_Expression_Estimation_Single_Sample_SE',
        'pipeline': 'RNASeqSingleSampleSE.xml',
        'option':
            'config_file_RSEM_RNA_SEQ_Single_Sample_Expression_Estimation',
        'args': fqse,
        'allow_groups': False,
        'testing_files': 'fqse',
    },
    'hrse_en': {
        'fullname': 'RNA Single End Engrafted tumor',
        'dir': 'Human/Xenome/RNA/RNA_Expression_Estimation_Single_Sample_SE',
        'pipeline': 'XenomeRnaSeqSingleSampleSE.xml',
        'option':
            'config_file_RSEM_RNA_SEQ_Single_Sample_Expression_Estimation',
        'args': fqse,
        'allow_groups': False,
        'testing_files': 'fqse',
        'help_prolog': """RNA Single End Engrafted tumor(hg38 version)""",
    },
    'hrpe_top': {
        'fullname':
            'Tophat RNA Paired End Patient only (variant calling too)',
        'dir': 'Human/RNA/TophatBased_RNA_Novel_isoform_detection_PE',
        'pipeline': 'Tophat_HumanRNASeqSingleSamplePE.xml',
        'option': 'config_file_RNA_SEQ_Single_Sample_Expression_Estimation',
        'args': fqpeg,
        'num_in_group': 2,
        'allow_groups': True,
        'testing_files': 'fqpeg',
        'help_prolog': """Tophat RNA Paired End Patient only (variant calling too) hg19""",
    },
    'hrrbs_pe': {
        'fullname':
            'Human RRBS Paired End by Bismark',
        'dir': 'Human/Bisulphite/Human_RRBS_PE_Bismark',
        'pipeline': 'Human_Bisulphite_RRBS_hg19.xml',
        'option': 'config_file_bisulphite_sequencing',
        'args': fqpeg,
        'num_in_group': 2,
        'allow_groups': True,
        'testing_files': 'fqpeg',
        'help_prolog': """ Human RRBS Paired End by Bismark (hg19)""",
    },
    'hrse_emase_en': {
        'fullname':
            'RNA-Seq Expression Estimation by EMASE for PDX Samples (Single end fastq file)',
        'dir': 'Human/Xenome/RNA/RNA_Expression_EMASE_SE',
        'pipeline': 'EMASE_PDX_RNA_SE.xml',
        'option': 'config_file_EMASE_param_PDX',
        'args': fqse,
        'allow_groups': False,
        'testing_files': 'fqse',
        'help_prolog': """RNA-Seq Expression Estimation using EMASE for PDX Samples (engrafted in Mice).
Please use single-end fastq as an input for the pipeline. Pipeline also accepts compressed formats(gz or bzip2)
This pipeline is for RNA Expression Estimation for PDX Samples(engrafted in Mice). 
Pipeline does the simultaneous alignment to human (hg38) and mouse (NOD_ShiLtJ based on mm10) transcriptome and expression estimation using EMASE.
The gtf version for Human and Mouse is  Homo_sapiens.GRCh38.84.chr_patch_hapl_scaff.gtf and Mus_musculus.GRCm38.82.chr.gtf,respectively.
The mouse transcripts and gene in the output file can be identified by "MUS" in the name.
""",
        'help_epilog': '''Example usage :
hrse_emase_en  -o config_file_EMASE_param_PDX  in_R1.fastq '''
    },
    'hrpe_fuse': {
        'fullname': 'RNA-seq fusion gene detection using Soapfuse',
        'dir': 'Human/Xenome/RNA/RNA_Soapfuse',
        'pipeline': 'XenomeRnaSeqSoapfuse.xml',
        'option': None,
        'args': fqpeg,
        'num_in_group': 2,
        'allow_groups': True,
        'testing_files': 'fqpeg',
        'help_prolog': """RNA-seq fusion gene detection using Soapfuse(hg38)""",
    },
    'hte_dfci_en': {
        'fullname':
            'Targeted exome for Paired end data Engrafted Tumor from DFCI',
        'dir': 'Human/Xenome/Exome/CCP',
        'pipeline': 'XenomeSingleSampleCCP.xml',
        'option': 'config_file_Xenome_DFCI_TE_Single_Sample',
        'args': fqpeg,
        'num_in_group': 2,
        'allow_groups': True,
        'testing_files': 'fqpeg',
        'help_prolog': """Targeted exome for Paired end data Engrafted Tumor from DFCI(hg38)""",

    },
    'hwes_baylor_en': {
        'fullname':
            'Whole exome for Paired end data Engrafted Tumor from Baylor',
        'dir': 'Human/Xenome/Exome/WholeExome',
        'pipeline': 'XenomeWholeExomeSingleSample.xml',
        'option': 'Baylor_config_file_Whole_Exome_Single_Sample',
        'args': fqpeg,
        'num_in_group': 2,
        'allow_groups': True,
        'testing_files': 'fqpeg',
        'help_prolog': """Whole exome for Paired end data Engrafted Tumor from Baylor(hg38)""",
    },
    'hmicr_16S_A': {
        'fullname':
            'Human Microbiome pipeline A for 16S data',
        'dir': 'Human/Microbiome/16SrRNA/1.0',
        'pipeline': 'pipeline-16S-A.xml',
        'option': 'config_file_16SPipelineA',
        'args': microbiome_dir,
        'num_in_group': 2,
        'allow_groups': False,
        'testing_files': 'microbiome_dir',
        'help_prolog': """Human Microbiome 16ssrRNA pipeline-A""",
        'help_epilog': '''Example usage :
          hmicr_16s_A   input_directory   output_directory  (command specifying path of input and output directories)'''

    },
    'mchip': {
        'fullname': 'ChIP-seq Paired End Experimental-Control',
        'dir': 'Mouse/ChIPseq/ChIPseqPairedEnd',
        'pipeline': 'experimentalControlChIP.xml',
        'option': None,
        'args': chipc,
        'allow_groups': False,
        'testing_files': 'chipc',
    },
    'hte_truseq_en': {
        'fullname': 'Targeted Sequencing Paired End data Trueseq panel Engrafted tumor',
        'dir': 'Human/Xenome/Exome/truseq_snpeff_exome',
        'pipeline': 'xenome_single_sample_exome_quadraploid.xml',
        'option': 'config_file_Xenome_truseq_snpeff_exome',
        'num_in_group': 2,
        'args': fqpe,
        'allow_groups': True,
        'testing_files': 'fqpe',
        'help_prolog': """Targeted Sequencing Paired End data Trueseq panel Engrafted tumor(hg38)""",
    },
    'hte_trusight_en': {
        'fullname': 'Targeted Sequencing Paired End data Trusight myeloid panel Engrafted tumor',
        'dir': 'Human/Xenome/Exome/truseq_snpeff_exome',
        'pipeline': 'xenome_single_sample_exome_quadraploid.xml',
        'option': 'config_file_Xenome_trusight_snpeff_exome',
        'num_in_group': 2,
        'args': fqpe,
        'allow_groups': True,
        'testing_files': 'fqpe',
        'help_prolog': """Targeted Sequencing Paired End data Trusight myeloid panel Engrafted tumor(hg38)""",
    },
    'hmirbase_exp': {
        'fullname': 'Human miRNA-precursor (miRBase release 21) based pipeline',
        'dir': 'Human/RNA/HmiRBase_Expression_SE',
        'pipeline': 'miRNA_HumanAnalysisFormMiRBase_SE.xml',
        'option': 'config.file',
        'args': fqse,
        'allow_groups': False,
        'testing_files': 'fqse',
        'help_prolog': """Human miRNA-precursor (miRBase release 21) based pipeline.This Pipeline required a unpaired fastq file along with a config file """,
    },
    'mwes_closest': {
        'fullname':
            'Whole Exome Paired End data with closest strain for BQSR',
        'dir': 'Mouse/Exome/WholeExomeWithClose',
        'pipeline': 'Mouse_WholeExomeWithClose.xml',
        'option': 'config_file_Whole_Exome_MMR',
        'args': fqclosest,
        'trailing_args': 2,
        'num_in_group': 2,
        'allow_groups': True,
        'testing_files': 'fqclosest',
    },
    'mwes': {
        'fullname':
            'Whole Exome Paired End data without closest strain',
        'dir': 'Mouse/Exome/WholeExomeNoClose',
        'pipeline': 'Mouse_WholeExome_NoClose.xml',
        'option': 'config_file_Whole_Exome_MMR',
        'args': fqpeg,
        'num_in_group': 2,
        'allow_groups': True,
        'testing_files': 'fqpeg',
    },
    'mwgs_closest': {
        'fullname':
            'Mouse Whole Genome Paired End data with closest strain '
            'for BQSR',
        'dir': 'Mouse/WholeGenome/WholeGenomeWithClose',
        'pipeline': 'Mouse_WholeGenome_WithClose.xml',
        'option': None,
        'args': fqclosest,
        'num_in_group': 2,
        'trailing_args': 2,
        'allow_groups': True,
        'testing_files': 'fqclosest',
    },
    'mwgs': {
        'fullname': 'Mouse Whole Genome Paired End data without closest strain',
        'dir': 'Mouse/WholeGenome/WholeGenomeNoClose',
        'pipeline': 'Mouse_WholeGenome_NoClose.xml',
        'option': None,
        'args': fqpeg,
        'num_in_group': 2,
        'allow_groups': True,
        'testing_files': 'fqpeg',
    },
    'mrpe': {
        'fullname': 'Mouse RNA Paired End',
        'dir': 'Mouse/RNA/RNA_Expression_Estimation_Single_Sample_PE',
        'pipeline': 'Mouse_RNASeqSingleSamplePE.xml',
        'option':
            'config_file_RSEM_RNA_SEQ_Single_Sample_Expression_Estimation',
        'args': fqpeg,
        'num_in_group': 2,
        'allow_groups': True,
        'testing_files': 'fqpeg',
    },
    'mrse': {
        'fullname': 'Mouse RNA Single End',
        'dir': 'Mouse/RNA/RNA_Expression_Estimation_Single_Sample_SE',
        'pipeline': 'Mouse_RNASeqSingleSampleSE.xml',
        'option':
            'config_file_RSEM_RNA_SEQ_Single_Sample_Expression_Estimation',
        'args': fqse,
        'allow_groups': False,
        'testing_files': 'fqse',
    },
    'mrpe_DOEMASE': {
        'fullname': 'Diversity Outbred Mouse RNAseq Paired End EMASE based Pipeline',
        'dir': 'Mouse/RNA/MRPE_DOEMASE',
        'pipeline': 'Pipeline_MRPE_DOEMASE.xml',
        'option': 'config_file_Pipeline_MRPE_DOEMASE',
        'args': fqpe,
        'allow_groups': False,
        'testing_files': 'fqpe',
    },
    'mrse_DOEMASE': {
        'fullname': 'Diversity Outbred Mouse RNAseq Single End EMASE based Pipeline',
        'dir': 'Mouse/RNA/MRSE_DOEMASE',
        'pipeline': 'Pipeline_MRSE_DOEMASE.xml',
        'option': 'config_file_Pipeline_MRSE_DOEMASE',
        'args': fqse,
        'allow_groups': False,
        'testing_files': 'fqse',
    },
    'mrpe_MREMASE': {
        'fullname': 'Mouse RNAseq Paired End Multi Reference EMASE',
        'dir': 'Mouse/RNA/MRPE_MREMASE',
        'pipeline': 'Pipeline_MRPE_MREMASE.xml',
        'option': 'config_file_Pipeline_MRPE_MREMASE',
        'args': fqpe,
        'allow_groups': False,
        'testing_files': 'fqpe',
    },
    'mrse_MREMASE': {
        'fullname': 'Mouse RNAseq Single End Multi Reference EMASE',
        'dir': 'Mouse/RNA/MRSE_MREMASE',
        'pipeline': 'Pipeline_MRSE_MREMASE.xml',
        'option': 'config_file_Pipeline_MRSE_MREMASE',
        'args': fqse,
        'allow_groups': False,
        'testing_files': 'fqse',
    },
    'mrpe_SREMASE': {
        'fullname': 'Mouse RNAseq Paired End Single Reference EMASE',
        'dir': 'Mouse/RNA/MRPE_SREMASE',
        'pipeline': 'Pipeline_MRPE_SREMASE.xml',
        'option': 'config_file_Pipeline_MRPE_SREMASE',
        'args': fqpe,
        'allow_groups': False,
        'testing_files': 'fqpe',
    },
    'mrse_SREMASE': {
        'fullname': 'Mouse RNAseq Single End Single Reference EMASE',
        'dir': 'Mouse/RNA/MRSE_SREMASE',
        'pipeline': 'Pipeline_MRSE_SREMASE.xml',
        'option': 'config_file_Pipeline_MRSE_SREMASE',
        'args': fqse,
        'allow_groups': False,
        'testing_files': 'fqse',
    },
    'mrse_DOGBRS': {
        'fullname': 'Diversity Outbred Mouse RNAseq Single End GBRS algorithm based Pipeline',
        'dir': 'Mouse/RNA/MRSE_DOGBRS',
        'pipeline': 'Pipeline_MRSE_DOGBRS.xml',
        'option': 'config_file_Pipeline_MRSE_DOGBRS',
        'args': fqse_gensex, 
        'allow_groups': False,
        'testing_files': 'fqse_gensex',
        'help_prolog': "This pipeline input is a single end fastq file followed by a string describing the generation and sex e.g., 'G32.F'"
    },
    'mrpe_DOGBRS': {
        'fullname': 'Diversity Outbred Mouse RNAseq Paired End GBRS algorithm based Pipeline',
        'dir': 'Mouse/RNA/MRPE_DOGBRS',
        'pipeline': 'Pipeline_MRPE_DOGBRS.xml',
        'option': 'config_file_Pipeline_MRPE_DOGBRS',
        'args': fqpe_gensex,
        'allow_groups': False,
        'testing_files': 'fqpe_gensex',
        'help_prolog': "This pipeline input is paired end fastq files followed by a string describing the generation and sex e.g., 'G32.F'"
    },
    'mrpe_top': {
        'fullname': 'Mouse Tophat RNA Paired End Only',
        'dir': 'Mouse/RNA/TophatBased_RNA_Novel_isoform_detection_variant_calling_PE',
        'pipeline': 'Tophat_MouseRNASeqSingleSamplePE.xml',
        'option': 'config_file_RNA_SEQ_Single_Sample_Expression_Estimation',
        'args': fqpeg,
        'num_in_group': 2,
        'allow_groups': True,
        'testing_files': 'fqpeg',
    },
    'mrrbs_pe': {
        'fullname':
            'Mouse RRBS Paired End by Bismark',
        'dir': 'Mouse/Bisulphite/Mouse_RRBS_PE_Bismark',
        'pipeline': 'Mouse_Bisulphite_RRBS_mm10.xml',
        'option': 'config_file_bisulphite_sequencing',
        'args': fqpe,
        'allow_groups': False,
        'testing_files': 'fqpe',
    },
    'hcnv': {
        'fullname': 'SNP Array CNV Single Sample',
        'dir': 'arrayBased/Human/CNV',
        'pipeline': 'single_sample_cnv_snparray.xml',
        'option': None,
        'args': cel_single,
        'allow_groups': False,
        'testing_files': 'cel_single',
        'help_prolog': """SNP Array CNV Single Sample(hg38)""",
    },
    'hcnv_tn': {
        'fullname': 'SNP Array CNV for Paired Tumor/Normal Samples',
        'dir': 'arrayBased/Human/CNV',
        'pipeline': 'paired_sample_cnv_snparray.xml',
        'option': None,
        'args': cel_paired,
        'allow_groups': False,
        'testing_files': 'cel_paired',
        'help_prolog': """SNP Array CNV for Paired Tumor/Normal Samples(hg38)""",
    },
    'mmirbase_exp': {
        'fullname': 'Mouse miRNA-precursor (miRBase release 21) based pipeline',
        'dir': 'Mouse/RNA/MmiRBase_Expression_SE',
        'pipeline': 'miRNA_MouseAnalysisFormMiRBase_SE.xml',
        'option': 'config.file',
        'args': fqse,
        'allow_groups': False,
        'testing_files': 'fqse',
        'help_prolog': """Mouse miRNA-precursor (miRBase release 21) based pipeline.This Pipeline required a unpaired fastq file along with a config file """,
    },
}
