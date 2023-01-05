This pipeline uses GBRS suite of tools for reconstructing genomes using RNA-Seq data 
from multiparent population and quantifying allele specific expression.The GBRS 
suite of tools is developed by Churchill lab at the Jackson Laboratory.

The pipeline uses 1 single end fastq.gz file as the input along with the 
generation and sex information. For example for Generation 23 and Females,
use the following command-

mrpe_DOGBRS reads.fastq.gz G23.F

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



## References:

Please refer to the following website for more information
http://churchill-lab.github.io/gbrs/