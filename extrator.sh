#!/bin/bash

snp="NA"
vcf="NA"
type="DOSAGE"
ind="NA"
ld=0.0
output="SNPsExtracted"


while getopts ":s:v:t:i:l:o:" opt; do 
    case $opt in 
	s) #SNP <FILE>
	    snp=$OPTARG
	    ;;
	v) #VCFs <PATH>
	    vcf=$OPTARG
	    ;;
	t) #Type <DOSAGE,GENOTYPE,LIKELIHOOD>
	    type=$OPTARG
	    ;;
	i) #individuals <FILE>
	    ind=$OPTARG
	    ;;
	l) #LD <FLOAT>
	    ld=$OPTARG
	    ;;
	o) #output <STRING>
	    output=$OPTARG
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG" >&2
	    exit 1
	    ;;
	:)
	    echo "Option -$OPTARG requires an argument. See README.md" >&2
	    exit 1
	    ;;
    esac
done

#If statement to check that not both cohort and individuals is used
## Define value to determine if cohort or individual is used

echo "Script called with the following arguments:"
echo "snp: $snp"
echo "vcf: $vcf"
echo "type: $type"
echo "indiduals: $ind"
echo "ld: $ld"
echo "output: $output"

if [ $ind = 'NA' ]; then
    echo 'no list of individuals specified. Extracting SNPs from all individuals in dataset'
fi

# Captures missing arguments:
if [ $snp = 'NA' ]; then
    echo 'no list of SNPs specified. Exiting..'
    exit 1
fi


#Setup phase
#currentExtract=`date | tr ' ' '.'`
currentExtract="Thu.Dec.15.15:34:27.CET.2016"
mkdir -p tmpGeno/imputed
mkdir -p tmpGeno/$currentExtract/extracted
mkdir -p tmpGeno/$currentExtract/splitted


if [ $vcf = 'NA' ]; then
    tmpvcf=`ls tmpGeno/imputed/`
    vcf=`echo "tmpGeno/imputed/$tmpvcf"`
    echo "no link to vcf given - using existing link at $vcf"
else
    if [ `echo "$vcf" | rev | cut -d'/' -f2 | rev` != "`ls tmpGeno/imputed`" ]; then
	ln -s $vcf tmpGeno/imputed/
	tmpvcf=`ls tmpGeno/imputed/`
	vcf=`echo "tmpGeno/imputed/$tmpvcf"`
    else
	echo "$vcf link already there"
	tmpvcf=`ls tmpGeno/imputed/`
	vcf=`echo "tmpGeno/imputed/$tmpvcf"`
    fi
fi

#splitting snplist into chromsomes, which is later looped through, for better I/O optimization (vcfs are not striped)
for chr in `awk '{print $2}' $snp | sort -n | uniq`
do
    awk -v "chr=$chr" '$2==chr {print $2"\t"$3"\t"$3}' $snp > tmpGeno/$currentExtract/splitted/SNPs.on.chr.$chr.list
done

if [ x = y ]; then
for snpchr in `ls tmpGeno/$currentExtract/splitted/SNPs.on.chr*.list`; do
    chr=`basename $snpchr | cut -d'.' -f4`
    if [ $ind = 'NA' ]; then
	echo "bcftools view -Oz -r $snpchr $vcf/$chr.vcf.gz -o tmpGeno/$currentExtract/extracted/from.chr.$chr.vcf.gz" 
    else
	echo "bcftools view -Oz -s $ind -r $snpchr $vcf/$chr.vcf.gz -o tmpGeno/$currentExtract/extracted/from.chr.$chr.vcf.gz" 	
    fi
done | ./submit_jobarray.py -m 18G -n logs.for.extract. #18G mem is based on size of chr2 for decode.
fi 
  
###check with snp-name?

#cat/collect everything into one vcf
subVCFs=`ls tmpGeno/$currentExtract/extracted/from.chr.*.vcf.gz | tr '\n' ' '`
echo "bcftools concat $subVCFs -Oz -o tmpGeno/$currentExtract/all.SNPs.extracted.vcf.gz" | ./submit_jobarray.py -w logs.for.extract. -n logs.for.concat. -m 10G

if [ $type = 'LIKELIHOOD' ]; then
    echo $type
elif [ $type = 'GENOTYPE' ]; then
    echo $type
else 
    echo $type
#    echo "bcftools annotate -R QUAL,FILTER,INFO,FORMAT/GT,FORMAT/ADS,FORMAT/GP tmpGeno/all.SNPs.extracted.vcf.gz" 
fi

#LD threshold

#Format into csv with R.
