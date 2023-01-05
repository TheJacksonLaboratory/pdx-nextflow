#! /usr/bin/env Rscript

library(ChIPpeakAnno)
library(EnsDb.Hsapiens.v79)
data("TSS.human.GRCh38")

#Input from chip_peak_anno.xml

args <- commandArgs(trailingOnly=TRUE)


chipanno <- function(args){
  macs <- read.csv(args[1], sep = "\t", col.names = c("chr","start","end","name","score"))
  #convert to GRanges, annotate with IDs. 
  macsOutput <- toGRanges(macs, format="MACS")
  macs.anno <- annotatePeakInBatch(macsOutput, AnnotationData=TSS.human.GRCh38, maxgap=args[3])
  #Preferred method: trim unused columns, annotate with gene symbol, entrez ID, and gene biotype.
  edb <- EnsDb.Hsapiens.v79
  macs.annoDF <- as.data.frame(macs.anno)
  macs.annoDF <- macs.annoDF[,!names(macs.annoDF) %in% "ensembl"]
  ensembl.anno <- genes(edb, filter = list(GeneIdFilter(macs.annoDF$feature)))
  ensembl.anno.DF <- as.data.frame(ensembl.anno)
  ensembl.anno.DF.cut <- ensembl.anno.DF[!names(ensembl.anno.DF) %in% c("seqnames","start","end","width", "strand", "seq_coord_system")]
  result <- merge(macs.annoDF, ensembl.anno.DF.cut, by.x="feature", by.y = "gene_id")
  ####Alternative method: ChIPpeakAnno-based annotation: Slightly more informative (contains full gene names), but fails to annotate the majority of predicted genes.
  #library("org.Hs.eg.db")
  #macs.annoL <- addGeneIDs(macs.anno,"org.Hs.eg.db",c("symbol", "genename","entrez_id"))
  #result <- as.data.frame(macs.annoL)
  #write.table(x=result, file=args[4], sep = "\t", quote = F, row.names = F)
  # writting annotation results with a comment
  con <- file(args[4], open="wt")
  writeLines(paste("#Peaks that do not correspond to canonical Ensembl-named loci (unlocalized and unplaced chr) are not included for annotation."), con)
  write.table(result, file=con, sep = "\t", quote = F, row.names = F)
  close(con)
}

chipanno(args)
