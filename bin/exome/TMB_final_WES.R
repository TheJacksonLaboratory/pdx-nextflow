#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

data <- read.table(args[1], header=T)

d<-data[!((data$chr== "chr3" & data$start=="195000000") | (data$chr == "chr11" & data$start=="1000000")  |  (data$chr == "chr19" & data$start == "8000000") |        (data$chr == "chr7" & data$start == "101000000") | (data$chr == "chr6" & data$start== "29000000")),]

print(sum(d$HM))
print(sum(d$length))

TMB<- (sum(d$HM)/ sum(d$length))*1000000
 
write.table(TMB, file=args[2], quote=F, row.names=F, col.names=F)


