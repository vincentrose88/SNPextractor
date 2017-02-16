#!/usr/bin/Rscript
args <- commandArgs(trailingOnly=T)

folder <- args[1]
outputName <- args[2]
ldcut <- args[3]

geno <- read.table(paste(folder,'genoFile.noHead',sep=''),h=F,as.is=T)
genoHeader <- read.table(paste(folder,'newHeader',sep=''),h=F,as.is=T)
colnames(geno) <- t(genoHeader)
#Writing info as numeric
geno$INFO <- as.numeric(substr(geno$INFO,6,nchar(geno$INFO)))

if(FALSE){
if(!is.na(ldcut)){
    ldfile <- read.table(paste(folder,'ld/plink.ld',sep=''),h=T,as.is=T)
    if(nrow(ldfile==0)){
        print(paste('No SNPs in LD, none removed due to ld-threshold of R2 <',ldcut))
    }else{
       ld <- ldfile[ldfile$R2>ldcut,c('SNP_A','SNP_B')]

       
       tmpmerge <- merge(geno[geno$ID %in% c(ld$SNP_A,ld$SNP_B),c('ID','INFO')],ld,by.x='ID',by.y='SNP_A')
#       colnames(tmpmerge)[colnames(tmpmerge)=='INFO'] <- 'INFO_A'
       ldmerge <- merge(geno[geno$ID %in% c(ld$SNP_A,ld$SNP_B),c('ID','INFO')],ld,by.x='ID',tmpmerge,by.y='SNP_B',all=T)
#       colnames(ldmerge)[colnames(ldmerge)=='INFO'] <- 'INFO_B'
       ldmerge <- ldmerge[,1:2]
    }
}
}

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

#Write out as csv's
write.csv(finalGeno,paste(outputName,'geno','csv',sep='.'),row.names=T,quote=F)
write.csv(finalInfo,paste(outputName,'info','csv',sep='.'),row.names=F,quote=F)
