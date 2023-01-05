This pipeline uses GBRS suite of tools for reconstructing genomes using RNA-Seq data 
from multiparent population and quantifying allele specific expression.The GBRS 
suite of tools is developed by Churchill lab at the Jackson Laboratory.

The pipeline uses paired-end fastq.gz files as the input along with the 
generation and sex information. For example for Generation 23 and Female,
use the following command-

mrpe_DOGBRS R1reads.fastq.gz R2reads.fastq.gz G23.F

It produces the following resulting files

* .gbrs.plotted.genome.pdf
* .gbrs.interpolated.genoprobs.npz
* .diploid.genes.alignment_counts
* .diploid.isoforms.alignment_counts
* .diploid.genes.expected_read_counts
* .diploid.isoforms.expected_read_counts
* .diploid.genes.tpm
* .diploid.isoforms.tpm
* .genotypes.tsv
* .genotypes.npz
* .genoprobs.npz
* .multiway.genes.alignment_counts
* .multiway.isoforms.alignment_counts
* .multiway.genes.expected_read_counts
* .multiway.genes.tpm
* .multiway.isoforms.expected_read_counts


In addition to the above files, this workflow also runs several QC steps:

1. **Untrimmed** reads are aligned against B6 genome with STAR, and mapping stats are collected.   
	>**NOTE** Alignment metrics are for information only. STAR alignment is not used in GBRS calculations.  
	
2. Picard `collectRNAseqMetrics` is run on the STAR alignment to calculate proportion ribosomal reads, and proportion of useable bases.  
 	>**NOTE** Metrics are for information only. 
 	
The QC metrics are aggregated into a single file: 

* \_summary\_stats.txt


## References:

Please refer to the following website for more information
http://churchill-lab.github.io/gbrs/