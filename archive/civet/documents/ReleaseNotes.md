# Release Notes: 

## 4.7.1

* Pipelines Added -
	* CTP Mutect2 based pipeline (pull request #32)
	
* Pipeline changes -
	* DO GBRS / EMASE (pull request #30)
		* Modified for use on Helix. 
		* Added RNA QC to PE\_GBRS pipeline
	* Bug fixes and minor changes for pipelines (pull request #27, pull request #29)

 * Pipeline removed -
	* Obsolete EMASE pipelines (pull request #31)

## 4.7.0
* Pipelines Added - 
	* DO EMASE (pull request #20)
		* MRPE_MREMASE 
		* MRPE_SREMASE
		* MRSE_MREMASE
		* MRSE_SREMASE

* Pipeline changes -
	* XML descriptions added (pull request #22)
		* Added metadata to pipelines for use with JODA.
		* Removed older e1\_fastq, remaining\_fastqs work-around.
		* Bug fixes for pipelines found during testing. 
	* Bug fixes for pipelines (pull request #19, pull request #21, pull request #23, pull request #24)

	
* Pipeline removed -
	* Mouse\_WholeGenome\_Singlecell

## 4.6.0
* Pipelines Added - 
	* DO EMASE/GBRS (pull request #11)
		* Diversity Outbred Mouse RNAseq Paired End EMASE based Pipeline
		* Diversity Outbred Mouse RNAseq Single End EMASE based Pipeline
		* Diversity Outbred Mouse RNAseq Paired End GBRS algorithm based Pipeline
		* Diversity Outbred Mouse RNAseq Single End GBRS algorithm based Pipeline
	* Metagenome Assembly (pull request #14)

* Pipeline Changes - 
	* General changes to pipelines (pull request #17):
		* Fixed AdjSNP problem (added extra step as XML tool)
		* CPP, WEX and TruSeq/TruSight pipelines now take a list of files as input.
		* Config files and their input are also fixed within all pipelines.
	* CNV pipeline (pull request #13)
		* Updated python script to get gender (get_model_gender.py)
		* Added CNV_seg plot
	* Conflicts in RNA analysis directory naming for DFCI and Baylor samples have been corrected. (pull request #15) 

* Pipelines removed -
	* Older version of GBRs pipelines (pull request #12)
		* Mouse 8-way alignment and perform Genotyping By RNA-Seq  (GBRS) Single end)
		* Mouse 8-way alignment and perform Genotyping By RNA-Seq  (GBRS) Paired end)
	