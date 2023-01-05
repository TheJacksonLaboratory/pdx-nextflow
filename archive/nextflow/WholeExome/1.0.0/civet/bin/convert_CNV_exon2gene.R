#! /usr/bin/env Rscript

##This is an R script for (1) recalibrating the exon-level p-values from CONTRA using conReg and (2) estimating gene level p-value from exon-level p-values using Fisher's method

args <- commandArgs(trailingOnly = TRUE)

rawinp=read.delim(file=args[1])
conreg_functions=args[2]

threshold_for_amplification=6     ##This is based on sensitivity data. Needs to be changed if there is an update to sensitivity analysis.

clinical=read.delim(file="/data/shared/cga_reference_data/clinical_targets_for_CNV.txt", header=T)

out_all_genes=paste(args[1],"_geneLevel_allGenes.txt",sep="")
out_sig_genes=paste(args[1],"_geneLevel_significantGenes.txt",sep="")
out_sig_genes_clinical=paste(args[1],"_geneLevel_significantGenes.clinicalTargets.cnv",sep="")

rawinp[is.na(rawinp$P.Value) == FALSE,] -> inp

#Recalibrating the exon-level p-values from CONTRA using conReg
source(conreg_functions)
p = inp$P.Value
obj = reg_pvalue(p)

cbind(obj$p_index,obj$p_input,obj$p_adj) -> obj1
obj1[order(obj1[,1]),]-> obj1s
obj1s[,3][which(obj1s[,3] < 0)] = 0 ##coerce negative p-values to 0

cbind(inp,obj1s[,3]) -> out
colnames(out) = c(colnames(inp), "P.Value.ConReg")

#Estimating gene level p-value from exon-level p-values using Fisher's method
get_gene_level_pval <- function(pvals) {   
k=-2*log2(prod(pvals))
return(pchisq(k,2*length(pvals),lower.tail=F))
}


cbind(aggregate(P.Value.ConReg ~ Gene.Sym+Chr, data=out, FUN=c), aggregate(Adjusted.Mean.of.LogRatio ~ Gene.Sym+Chr, data=out, FUN=mean)) -> agg

agg_pvalue_adjusted =matrix(nrow=nrow(agg),ncol=5)
for (row in 1:nrow(agg)){
agg_pvalue_adjusted[row,1]=as.character(agg[row,]$Gene.Sym)
agg_pvalue_adjusted[row,2]=as.character(agg[row,]$Chr)
agg_pvalue_adjusted[row,3]=as.character(agg[row,]$Adjusted.Mean.of.LogRatio)
agg_pvalue_adjusted[row,4]=as.character(2*2^agg[row,]$Adjusted.Mean.of.LogRatio)
agg_pvalue_adjusted[row,5]=get_gene_level_pval(unlist(agg[row,]$P.Value.ConReg))
}

colnames(agg_pvalue_adjusted)=c("Gene", "Chr", "Adjusted.Mean.of.LogRatio", "CN.State", "P.Value")
agg_pvalue_adjusted=data.frame(agg_pvalue_adjusted)
round(as.numeric(as.character(agg_pvalue_adjusted$CN.State)), 1) -> agg_pvalue_adjusted$CN.State
round(as.numeric(as.character(agg_pvalue_adjusted$Adjusted.Mean.of.LogRatio)), 2) -> agg_pvalue_adjusted$Adjusted.Mean.of.LogRatio
signif(as.numeric(as.character(agg_pvalue_adjusted$P.Value)), 3) -> agg_pvalue_adjusted$P.Value
agg_pvalue_adjusted$CN.Call="normal"
agg_pvalue_adjusted[round(agg_pvalue_adjusted$CN.State) >= threshold_for_amplification,"CN.Call"]="amplified"
agg_pvalue_adjusted$LOH="-"
colnames(agg_pvalue_adjusted)=c("Gene", "Chr", "Adjusted.Mean.of.LogRatio", "CN State", "P.Value", "CN Call", "LOH")
agg_pvalue_adjusted[agg_pvalue_adjusted$P.Value <= 0.05 & abs(agg_pvalue_adjusted$Adjusted.Mean.of.LogRatio) > 0.3,] -> agg_pvalue_adjusted_sig_genes

merge(agg_pvalue_adjusted_sig_genes, clinical, by.x="Gene", by.y="Gene") -> agg_pvalue_adjusted_sig_genes_clinical

write.table(agg_pvalue_adjusted,file=out_all_genes,append=F,quote=F,row.names=F,col.names=T,sep="\t")
write.table(agg_pvalue_adjusted_sig_genes,file=out_sig_genes,append=F,quote=F,row.names=F,col.names=T,sep="\t")
write.table(agg_pvalue_adjusted_sig_genes_clinical,file=out_sig_genes_clinical,append=F,quote=F,row.names=F,col.names=T,sep="\t")