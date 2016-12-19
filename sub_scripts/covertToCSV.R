#!/usr/bin/Rscript
args <- commandArgs(trailingOnly=T)

folder <- args[1]
outputName <- args[2]

geno <- read.table(paste(folder,'genoFile.noHead',sep=''),h=F,as.is=T)
genoHeader <- read.table(paste(folder,'newHeader',sep=''),h=F,as.is=T)
colnames(geno) <- t(genoHeader)

finalGeno <- geno[,-which(colnames(geno) %in% c("QUAL", "FILTER", "INFO", "FORMAT","REF","ALT","POS","CHROM"))]
colnames(finalGeno)[colnames(finalGeno)=="ID"] <- "SNP"
finalGeno <- t(finalGeno)
colnames(finalGeno) <- finalGeno['SNP',]
finalGeno <- finalGeno[-which(rownames(finalGeno)=='SNP'),]

finalInfo <- geno[,which(colnames(geno) %in% c("ID","REF","ALT","POS","CHROM"))]
colnames(finalInfo)[colnames(finalInfo)=="ID"] <- "SNP"

finalInfo <- finalInfo[,c('SNP','CHROM','POS','REF','ALT')]
colnames(finalInfo)[4:5] <- c("REF(0)","ALT(2)")

write.csv(finalGeno,paste(outputName,'geno','csv',sep='.'),row.names=T,quote=F)
write.csv(finalInfo,paste(outputName,'info','csv',sep='.'),row.names=F,quote=F)
