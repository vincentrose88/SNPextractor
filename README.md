# SNPextractor
Simple extractor script to extract list of SNPs from HRC VCFs for CBMR

##Requirement
 * Software/build upon: 
  * bcftools 1.3.1
  * R 3.0.1
  * python 2.7.3
Data
 * Imputed dataset on HRC panel (2016), in .vcf format, split into chromosomes
 * List of SNPs with chromosome and position on build 37 (same as HRC panel). One SNP pr. line in the following order: SNPID chromosome positions
 * List of indidviduals as particids. One individual pr. line.

##How to install
clone this repository with:
`git clone git@github.com:vincentrose88/SNPextractor.git`

which will create the folder `SNPextractor` from where you will run the program

##How to use

Run the extractor.sh script with the following flags and arguments:
 * `-s` (**S**NP): File containing the list of SNPs that needs to extracted with the following requirement:
  * Each line has variant name, chromosome and position, seperated by `\t` (tab) (see example file). 
  * **NB: Must be on same build as imputed data (commonly b37, but do check)**. 
 * `-v` (**V**CF): Path to vcf containing the HRC imputed data. 
  * For example: `/emc/cbmr/data/imputed/decode-Nov09-sanger-hwe10e5/`
 * `-t` (**t**ype): Specify which type of data for each SNP is needed. **Only one is allowed**. 
  * Options are: DOSAGE, GENOTYPE or LIKELIHOOD. Default: DOSAGE
  * `-i` (**i**ndividuals): File containing a list of individuals you want extracted. One individuals particid on each line, using "\n".
 * Optional parameters
  * `-l` **L**D cut-off (for R^2) for excluding variants which are in LD with eachothers in the dataset. Only the variant with the highest imputation quality (info score) will be kept.
  * `-o` **o**utput filename. The suffix `.csv` will be added regardless. Default is: SNPsExtracted.csv
  * `-d` **d**ate flag. Use if you want to have timestamp-names on the temporary work-directory. Mandatory when using the extracting script from the same root folder to avoid overwriting.
  * `-m` **m**emory requested post chromosomes extraction. How much memory do you want to request from the cluters grid engine after extracting from chromosomes? Default: 2G
## Example
Extract dosage information for SNPs not in LD (R^2 < 0.5) for all inter99 individuals from Decode dataset:

`./extractor.sh -s SNPlist.example -v /emc/cbmr/data/imputed/decode-Nov09-sanger-hwe10e5/ -t DOSAGE -c 1000 -l 0.5 -o inter99.decode.extracted`

