#!/bin/bash

snp="NA"
vcf="NA"
type="DOSAGE"
ind="NA"
ld="NA"
output="SNPsExtracted"
useDate=true
memory="2G"
cohort="NA"
multipleCohort="NA"
multipleCohortUsed=false
sanger=false
grs=false

while getopts ":s:v:t:i:l:o:m:c:dag" opt; do 
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
	m) #memory <STRING>
	    memory=$OPTARG
	    ;;
	c) #cohort <integer>
	    cohort=$OPTARG
	    ;;
	d) #date flag
	    useDate=false
	    ;;
	a) #sAnger server
	    sanger=true
	    ;;
	g) #grs output?
	    grs=true
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

#Readies the check for the final files

isTheOutputFileThere="$output.geno.csv"
isTheOutputInfoThere="$output.info.csv"

#Special GRS case
if $grs; then
    type="DOSAGE"
    ld="NA"
    isTheOutputFileThere="genoFile.noHead"
    isTheOutputInfoThere="newHeader"
fi

echo "Script called with the following arguments:"
echo "snp: $snp"
echo "vcf: $vcf"
echo "type: $type"
echo "indiduals: $ind"
echo "cohort studyid: $cohort"
echo "ld: $ld"
echo "output: $output"
echo "Use timestamp: $useDate"
echo "memory requested (post extracting from chromosomes): $memory"
echo "sanger-flag set to: $sanger"
echo "GRS-flag set to: $grs"
echo "-------------------------------------------------------------------------------------"

#Check if multiple cohorts is given
if [[ $cohort != 'NA' ]] && [[ $cohort == *","* ]]; then
    multipleCohort=`echo $cohort | sed 's/^/$2==/g' | sed 's/,/ || $2==/g' `
    multipleCohortUsed=true
fi


# Captures missing arguments:
if [[ $snp = 'NA' ]]; then
    echo 'no list of SNPs specified. Exiting..'
    exit 1
fi


#Setup phase
if $useDate; then
    currentExtract=`date | tr ' ' '.'`
else
    currentExtract='tmpWorkdir'
fi
mkdir -p tmpGeno/imputed
mkdir -p tmpGeno/$currentExtract/extracted
mkdir -p tmpGeno/$currentExtract/splitted

#Clearing out old output files
if [ -f "$isTheOutputFileThere" ]; then
    echo "Removing old outputfile: $output.csv"
    rm $isTheOutputFileThere
    rm $isTheOutputInfoThere
fi 

if [ $vcf = 'NA' ]; then
    tmpvcf=`ls tmpGeno/imputed/ | head -1`
    vcf=`echo "tmpGeno/imputed/$tmpvcf"`
    echo "no link to vcf given - using first existing link at $vcf"
else
    if [ `echo "$vcf" | rev | cut -d'/' -f2 | rev` != "`ls tmpGeno/imputed/ | head -1`" ]; then
	ln -s $vcf tmpGeno/imputed/
	tmpvcf=`ls tmpGeno/imputed/ | head -1`
	vcf=`echo "tmpGeno/imputed/$tmpvcf"`
    else
	echo "$vcf link already there"
	tmpvcf=`ls tmpGeno/imputed/ | head -1`
	vcf=`echo "tmpGeno/imputed/$tmpvcf"`
    fi
fi

if $sanger; then
    if [ "x$(ls $vcf* | grep dose)" != "x" ]; then
	echo "Sanger flag is off, but vcf-names are from Sanger. Have you forgotten the sanger-flag (-a)?"
	exit 3
    fi
else
    if [ "x$(ls $vcf* | grep dose)" == "x" ]; then
	echo "Sanger flag is on, but vcf-names are from Michigan. Have you the sanger-flag (-a) in your command?"
	exit 4
    fi
fi


#Use cohort studyid if there is no individual list.
if [[ $ind = 'NA' ]] && [[ $cohort = 'NA' ]]; then
    echo 'no list of individuals specified and no cohort study id given. Extracting SNPs from all individuals in dataset'
elif [[ $ind = 'NA' ]]; then 
    echo "using cohort study id(s): $cohort to extract individuals"
    #getting header from vcf to only get genotyped individuals - hardcorded to extract from smallest chromosome (22)
    bcftools view -h $vcf/22.vcf.gz | grep '#CHROM' | cut -f10- | tr '\t' '\n' > tmpGeno/$currentExtract/ind.in.vcf
    if $multipleCohortUsed; then
	cohortOut=`echo $cohort | sed 's/,/./g'`
	grep -wFf <(awk "$multipleCohort" ./sub_scripts/particid.studyid.list | cut -d' ' -f1) tmpGeno/$currentExtract/ind.in.vcf > studyid.$cohortOut.in.vcf.individuals.list
	ind="studyid.$cohortOut.in.vcf.individuals.list"
    else
	grep -wFf <(awk -v "cohort=$cohort" '$2==cohort {print $1}' ./sub_scripts/particid.studyid.list) tmpGeno/$currentExtract/ind.in.vcf > studyid.$cohort.in.vcf.individuals.list
	ind="studyid.$cohort.in.vcf.individuals.list"
    fi
fi



#splitting snplist into chromsomes, which is later looped through, for better I/O optimization (vcfs are not striped)
for chr in `awk '{print $2}' $snp | sort -n | uniq`
do
    awk -v "chr=$chr" '$2==chr {print $2"\t"$3"\t"$3}' $snp > tmpGeno/$currentExtract/splitted/SNPs.on.chr.$chr.list
done

for chr in `awk '{print $2}' $snp | sort -n | uniq`
do
    snpchr="tmpGeno/$currentExtract/splitted/SNPs.on.chr.$chr.list"
    if $sanger; then
	if [ $ind = 'NA' ]; then
	    echo "bcftools view -Oz -R $snpchr $vcf/$chr.vcf.gz -o tmpGeno/$currentExtract/extracted/from.chr.$chr.vcf.gz" 
	else
	    echo "bcftools view -Oz -S $ind -R $snpchr $vcf/$chr.vcf.gz -o tmpGeno/$currentExtract/extracted/from.chr.$chr.vcf.gz" 	
	fi
    else
	if [ $ind = 'NA' ]; then
	    echo "bcftools view -Oz -R $snpchr $vcf/chr$chr.dose.vcf.gz -o tmpGeno/$currentExtract/extracted/from.chr.$chr.vcf.gz" 
	else
	    echo "bcftools view -Oz -S $ind -R $snpchr $vcf/chr$chr.dose.vcf.gz -o tmpGeno/$currentExtract/extracted/from.chr.$chr.vcf.gz" 	
	fi
    fi
done | ./sub_scripts/submit_jobarray.py -m 18G -n extract. #18G mem is based on size of chr2 for decode.

#New tactic - pre-emptive assign subVCFs before extraction is done from input file - and make the concat wait for the extract script

tmpsubVCFs="$(for splitName in tmpGeno/$currentExtract/splitted/SNPs.on.chr*.list
do
     chr=`basename $splitName | cut -d'.' -f 4`
     echo "tmpGeno/$currentExtract/extracted/from.chr.$chr.vcf.gz"
done)"

subVCFs=`echo $tmpsubVCFs | tr '\n' ' '`

#Let the fileserver catch up:
#echo "sleep 30" | ./sub_scripts/submit_jobarray.py -n catchUp. -w extract.

#cat/collect everything into one vcf

#Concating
echo "bcftools concat $subVCFs -Oz -o tmpGeno/$currentExtract/all.SNPs.extracted.vcf.gz" | ./sub_scripts/submit_jobarray.py -n concat. -m $memory -w extract. #catchUp.

#Getting the header (for later use)
echo "bcftools view -h tmpGeno/$currentExtract/all.SNPs.extracted.vcf.gz -Ov -o tmpGeno/$currentExtract/header.vcf" | ./sub_scripts/submit_jobarray.py -n header. -w concat. -m $memory
echo "./sub_scripts/vcfToR.header.sh tmpGeno/$currentExtract"  | ./sub_scripts/submit_jobarray.py -w header. -n finalHeader. -m $memory

if $sanger; then
    if [ $type = 'LIKELIHOOD' ]; then
	echo "bcftools annotate -Oz -x QUAL,FILTER,INFO/RefPanelAF,INFO/AC,INFO/AN,FORMAT/GT,FORMAT/ADS,FORMAT/DS tmpGeno/$currentExtract/all.SNPs.extracted.vcf.gz -o tmpGeno/$currentExtract/all.SNPs.formatted.vcf.gz" | ./sub_scripts/submit_jobarray.py -m $memory -n format. -w concat.
    elif [ $type = 'GENOTYPE' ]; then
	echo "bcftools annotate -Oz -x QUAL,FILTER,INFO/RefPanelAF,INFO/AC,INFO/AN,FORMAT/DS,FORMAT/ADS,FORMAT/GP tmpGeno/$currentExtract/all.SNPs.extracted.vcf.gz -o tmpGeno/$currentExtract/all.SNPs.formatted.vcf.gz" | ./sub_scripts/submit_jobarray.py -m $memory -n format. -w concat.
    else 
	echo "bcftools annotate -Oz -x QUAL,FILTER,INFO/RefPanelAF,INFO/AC,INFO/AN,FORMAT/GT,FORMAT/ADS,FORMAT/GP tmpGeno/$currentExtract/all.SNPs.extracted.vcf.gz -o tmpGeno/$currentExtract/all.SNPs.formatted.vcf.gz" | ./sub_scripts/submit_jobarray.py -m $memory -n format. -w concat.
    fi
else
    if [ $type = 'LIKELIHOOD' ]; then
	echo "bcftools annotate -Oz -x QUAL,FILTER,INFO/AF,INFO/MAF,INFO/R2,INFO/ER2,FORMAT/GT,FORMAT/DS tmpGeno/$currentExtract/all.SNPs.extracted.vcf.gz -o tmpGeno/$currentExtract/all.SNPs.formatted.vcf.gz" | ./sub_scripts/submit_jobarray.py -m $memory -n format. -w concat.
    elif [ $type = 'GENOTYPE' ]; then
	echo "bcftools annotate -Oz -x QUAL,FILTER,INFO/AF,INFO/MAF,INFO/R2,INFO/ER2,FORMAT/DS,FORMAT/GP tmpGeno/$currentExtract/all.SNPs.extracted.vcf.gz -o tmpGeno/$currentExtract/all.SNPs.formatted.vcf.gz" | ./sub_scripts/submit_jobarray.py -m $memory -n format. -w concat.
    else 
	echo "bcftools annotate -Oz -x QUAL,FILTER,INFO/AF,INFO/MAF,INFO/R2,INFO/ER2,FORMAT/GT,FORMAT/GP tmpGeno/$currentExtract/all.SNPs.extracted.vcf.gz -o tmpGeno/$currentExtract/all.SNPs.formatted.vcf.gz" | ./sub_scripts/submit_jobarray.py -m $memory -n format. -w concat.
    fi
fi
echo "bcftools view -Ov -H  tmpGeno/$currentExtract/all.SNPs.formatted.vcf.gz -o tmpGeno/$currentExtract/genoFile.noHead"  | ./sub_scripts/submit_jobarray.py -m $memory -n noheader. -w format.

#LD threshold
if [[ $ld != 'NA' ]]; then
    mkdir -p tmpGeno/$currentExtract/ld
    echo "plink19 --vcf tmpGeno/$currentExtract/all.SNPs.extracted.vcf.gz --r2 --out tmpGeno/$currentExtract/ld/plink" | ./sub_scripts/submit_jobarray.py -m $memory -n ld. -w concat.
    ldFile="tmpGeno/$currentExtract/ld/plink.ld"   
#Format into csv with R.
    echo "./sub_scripts/covertToCSV.R tmpGeno/$currentExtract/ $output $ld" | ./sub_scripts/submit_jobarray.py -n convertToFinal. -w noheader. -m $memory
else
    if [ "$grs" = false ] ; then
	echo "./sub_scripts/covertToCSV.R tmpGeno/$currentExtract/ $output" | ./sub_scripts/submit_jobarray.py -n convertToFinal. -w noheader. -m $memory
    else
	echo "./sub_scripts/convertToGRS.sh tmpGeno/$currentExtract/" | ./sub_scripts/submit_jobarray.py -w noheader. -n convertToFinal. -m $memory
    fi  
fi

##Waiting for qstat scripts:
while ! [ -f $isTheOutputFileThere ]
do
    ls > /dev/null
    if [ -f $isTheOutputFileThere ]; then
	break
    fi
    echo "-------------------------------------------------------------------------------------"
    echo '(Still) waiting for the cluster-submitted scripts to finish:'
    qstat | awk 'NR!=2 {print $1,$3,$5}'
    echo ''
    echo 'Or waiting for results file to be written'
    sleep 10
done

##Clean up and finish
echo "-------------------------------------------------------------------------------------"
mkdir -p logs
for i in e command o exe
do
    if [[ $ld != 'NA' ]]; 
    then
	for j in extract concat header finalHeader format convertToFinal ld noheader
	do
	    mv $j.*.$i logs/
	done
    else
	for j in extract concat header finalHeader format convertToFinal noheader
	do
	    mv $j.*.$i logs/
	done
   fi
done

echo "Extraction done. See logs/ for logs for each step."
if $grs; then
    echo "newHeader and genoFile.noHead is ready for GRS"
else
    nrOfSNPs=$(expr `wc -l $output.info.csv | cut -d' ' -f1` - 1)
    nrOfIndividuals=$(expr `wc -l $output.geno.csv | cut -d' ' -f1` - 1)
    echo "Your extraction of $nrOfSNPs SNPs for $nrOfIndividuals individuals is recorded in $output.geno.csv and $output.info.csv" 
fi
echo "-------------------------------------------------------------------------------------"
exit 0
