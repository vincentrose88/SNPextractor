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

#Getting the final format for the genofile
finalGeno <- geno[,-which(colnames(geno) %in% c("QUAL", "FILTER", "INFO", "FORMAT","REF","ALT","POS","CHROM"))]
colnames(finalGeno)[colnames(finalGeno)=="ID"] <- "SNP"
finalGeno <- t(finalGeno)
genoColnames <- finalGeno['SNP',]

if(ncol(finalGeno)>1){
    finalGeno <- finalGeno[-which(rownames(finalGeno)=='SNP'),]
}else{
    finalGeno <- as.data.frame(finalGeno[-which(rownames(finalGeno)=='SNP'),],stringsAsFactors=F)
}
colnames(finalGeno) <- genoColnames

#Info file
finalInfo <- geno[,which(colnames(geno) %in% c("ID","REF","ALT","POS","CHROM","INFO"))]
colnames(finalInfo)[colnames(finalInfo)=="ID"] <- "SNP"
finalInfo <- finalInfo[,c('SNP','CHROM','POS','REF','ALT','INFO')]
colnames(finalInfo)[4:5] <- c("REF(0)","ALT(2)")

#Adding users SNPnames from SNPlist file:
SNPnames <- read.table(userSNPsFile,h=F,as.is=T)
colnames(SNPnames) <- c('ID','CHROM','POS')

info <- merge(finalInfo,SNPnames,by=c('CHROM','POS')


#Write out as csv's
write.csv(finalGeno,paste(outputName,'geno','csv',sep='.'),row.names=T,quote=F)
write.csv(info,paste(outputName,'info','csv',sep='.'),row.names=F,quote=F)
