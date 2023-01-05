## Mouse Whole Exome Analysis Pipeline:

This pipeline maps bisulfite treated sequencing reads to the mm10 mouse reference genome, and perform methylation calls using the tool [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/).


### Workflow steps

#### Step1 Quality Trimming

Reads from FASTQ files are quality trimmed using a python script wrapped around **Trim Galore**

#### Step 2: Alignment and methylation call reporting (Bismark)

Filtered reads are aligned using Bismark with `--bowtie2` specified. 

The process is described in depth in the Bismark documentation: 

>Sequence reads are first transformed into fully bisulfite-converted forward (C->T) and reverse read (G->A conversion of the forward strand) versions, before they are aligned to similarly converted versions of the genome (also C->T and G->A converted). Sequence reads that produce a unique best alignment from the four alignment processes against the bisulfite genomes (which are running in parallel) are then compared to the normal genomic sequence and the methylation state of all cytosine positions in the read is inferred. A read is considered to align uniquely if an alignment has a unique best alignment score (as reported by the AS:i field). If a read produces several alignments with the same number of mismatches or with the same alignment score (AS:i field), a read (or a read-pair) is discarded altogether.

