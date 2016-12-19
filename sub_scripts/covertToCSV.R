#!/usr/bin/Rscript
args <- commandArgs(trailingOnly=T)

folder <- args[1]
outputName <- args[2]

geno <- read.table(paste(folder,'genoFile.noHead',sep=''),h=F,as.is=T)
genoHeader <- read.table(paste(folder,'newHeader',sep=''),h=F,as.is=T)
colnames(geno) <- t(genoHeader)

final <- geno[,-which(colnames(geno) %in% c("QUAL", "FILTER", "INFO", "FORMAT"))]

write.csv(final,paste(outputName,'csv',sep='.'),row.names=F,quote=F)
