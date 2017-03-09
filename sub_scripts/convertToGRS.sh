cp $1/all.SNPs.extracted.vcf.gz all.SNPs.vcf
sed -r 's/([0-9]+)\t([0-9]+)\t\./\1\t\2\tchr\1\:\2/g' $1/genoFile.noHead > genoFile.noHead
mv $1/newHeader .
