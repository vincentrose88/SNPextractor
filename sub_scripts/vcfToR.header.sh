tail -1 $1/header.vcf | sed 's/#//g' | sed 's/\t/"\t"/g' | sed 's/$/"/g' | sed 's/^/"/g' > $1/newHeader
