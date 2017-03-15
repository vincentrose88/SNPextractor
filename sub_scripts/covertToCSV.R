#!/usr/bin/Rscript
args <- commandArgs(trailingOnly=T)

folder <- args[1]
outputName <- args[2]
userSNPsFile <- args[3]

geno <- read.table(paste(folder,'genoFile.noHead',sep=''),h=F,as.is=T)
genoHeader <- read.table(paste(folder,'newHeader',sep=''),h=F,as.is=T)
colnames(geno) <- t(genoHeader)
#Writing info as numeric
geno$INFO <- as.numeric(substr(geno$INFO,6,nchar(geno$INFO)))

#Info file
info <- geno[,which(colnames(geno) %in% c("ID","REF","ALT","POS","CHROM","INFO"))]

#Adding users SNPnames from SNPlist file:
SNPnames <- read.table(userSNPsFile,h=F,as.is=T)
colnames(SNPnames) <- c('SNP','CHROM','POS')

finalInfo <- merge(info,SNPnames,by=c('CHROM','POS'),all.x=T)

#Ordering
finalInfo <- finalInfo[,c('ID','SNP','CHROM','POS','REF','ALT','INFO')]
colnames(finalInfo)[5:6] <- c("REF(0)","ALT(2)")



#Getting the final format for the genofile
finalGeno <- geno[,-which(colnames(geno) %in% c("QUAL", "FILTER", "INFO", "FORMAT","REF","ALT","POS","CHROM"))]
finalGeno <- t(finalGeno)
genoColnames <- finalGeno['ID',]

if(ncol(finalGeno)>1){
    finalGeno <- finalGeno[-which(rownames(finalGeno)=='ID'),]
}else{
    finalGeno <- as.data.frame(finalGeno[-which(rownames(finalGeno)=='ID'),],stringsAsFactors=F)
}
colnames(finalGeno) <- genoColnames
colnames(finalGeno)[1] <- 'particid'

#Write out as csv
write.csv(finalGeno,paste(outputName,'geno','csv',sep='.'),row.names=T,quote=F)
write.csv(finalInfo,paste(outputName,'info','csv',sep='.'),row.names=F,quote=F)
