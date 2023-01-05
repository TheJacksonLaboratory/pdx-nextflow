## Mouse Whole Exome Analysis Pipeline:

This Whole Exome Sequencing (WES) pipeline identifies variants from a mouse whole exome experiment by primarily using the [Broad Institute's](https://software.broadinstitute.org/gatk/best-practices/) best-practices workflow for alignment and GATK-3 HaplotypeCaller for variant calling. 


### Reference Files and Workflow Details

Required reference input files:

1. Target region files (BED and INTERVALS\_LIST format)
The BED file should correspond to the data being processed. Chromosome naming in the BED file should match the reference FASTA used to process the data (if using the default input FASTA file, please ensure that chromosome names begin with 'chr'). If a corresponding INTERVALS\_LIST file (as used by Picard toolkit) is not available, it can easily be generated using the **GATK BedToIntervalList** tool.

### Workflow steps and notable parameters

#### Step1 **Alignment and Target Coverage** 

Reads from FASTQ files trimmed (**Python Script**), aligned (**BWA**, alt-aware), sorted (**Picard SortSam**), and prepared for variant calling (**Picard MarkDuplicates**, **GATK BaseRecalibrator**, **GATK ApplyBQSR**). QC is also performed on the processed BAM files, and **Picard CalculateHsMetrics**. QC reports are aggregated with the **QC Integrate** script.

#### Step 2: Indexing BAM files (**Samtools Index BAM**)

#### Step 3: Variant calling and filtration (**GATK 3 HaplotypeCaller** and **VariantFiltration**)

#### Step 4: Annotation 

Basic annotation is done with **SnpEff** (using mm10 database). 
**dbNSFP** is used to add additional annotations from the [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP) database. Highest impact snpEff annotations are reported. 