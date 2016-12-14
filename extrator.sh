#!/bin/bash

snp="NA"
vcf="NA"
type="DOSAGE"
cohort="NA"
ind="NA"
ld=0.0
output="SNPsExtracted"


while getopts ":s:v:t:c:i:l:o:" opt; do 
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
	c) #cohort <prefix>
	    cohort=$OPTARG
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
echo "cohort: $cohort"
echo "indiduals: $ind"
echo "ld: $ld"
echo "output: $output"

mkdir -p tmpGeno
#If vcfs not already present:
ln -s $vcf tmpGeno/

#While read snp list
##bcftools position > tmp-file
###check with snp-name?

#cat/collect everything into one vcf
##bcftools <type> and <individuals/cohort> to only keep relevant info

#Format into csv with R.
