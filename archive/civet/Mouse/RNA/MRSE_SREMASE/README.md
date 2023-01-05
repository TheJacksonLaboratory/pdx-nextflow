This pipeline uses Expectation-Maximization (EMASE) algorithm for Allele
Specific Expression of mice with single reference genome. The EMASE algorithm
is developed by Churchill lab at the Jackson Laboratory.
The pipeline uses single end fastq.gz file as the input and produces the
following resulting files

* .isoforms.expected_read_counts
* .isoforms.tpm
* .genes.expected_read_counts
* .genes.tpm

Please refer to the following website for more information
https://emase.readthedocs.io/en/latest/readme.html

## References:
Narayanan Raghupathy, Kwangbom Choi, Matthew J. Vincent, Glen L. Beane, 
Keith Sheppard, Steven C. Munger, Ron Korstanje, 
Fernando Pardo-Manual de Villena, and Gary A. Churchill. 
Hierarchical Analysis of RNA-Seq Reads Improves the Accuracy of 
Allele-specific Expression Bioinformatics,bty078, 2018


Steven C. Munger, Narayanan Raghupathy, Kwangbom Choi, Allen K. Simons, 
Daniel M. Gatti, Douglas A. Hinerfeld, Karen L. Svenson, Mark P. Keller, 
Alan D. Attie, Matthew A. Hibbs, Joel H. Graber, Elissa J. Chesler 
and Gary A. Churchill.RNA-Seq Alignment to Individualized Genomes 
Improves Transcript Abundance Estimates in Multiparent Populations 
Genetics. 2014 Sep;198(1):59-73. doi: 10.1534/genetics.114.165886

Christopher L. Baker, Shimpei Kajita, Michael Walker, Ruth L. Saxl, 
Narayanan Raghupathy, Kwangbom Choi, Petko M. Petkov, Kenneth Paigen 
PRDM9 Drives Evolutionary Erosion of Hotspots in Mus musculus through 
Haplotype-Specific Initiation of Meiotic Recombination
PLOS Genetics: 2015;doi/10.1371/journal.pgen.1004916