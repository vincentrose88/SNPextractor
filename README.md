# SNPextractor
Simple extractor script to extract list of SNPs from HRC VCFs for CBMR

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
 * Either of the following but not both.
  * `-c` (**c**ohort**): Prefix for a specific cohort. Extract all individuals which matches this prefix. For example `57x` matches the CIMT cohort
  * `-i` (**i**ndividuals): File containing a list of individuals you want extracted. One individuals particid on each line, using "\n".
 * Optional parameters
  * `-l` **L**D cut-off (for R^2) for excluding variants which are in LD with eachothers in the dataset. Only the variant with the highest imputation quality (info score) will be kept.
  * `-o` **o**utput filename. The suffix `.csv` will be added regardless. Default is: SNPsExtracted.csv

## Example
Extract dosage information for SNPs not in LD (R^2 < 0.5) for all inter99 individuals from Decode dataset:

`./extractor.sh -s SNPlist.example -v /emc/cbmr/data/imputed/decode-Nov09-sanger-hwe10e5/ -t DOSAGE -c 1000 -l 0.5 -o inter99.decode.extracted`

