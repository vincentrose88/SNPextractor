# SNPextractor
Simple extractor script to extract list of SNPs from HRC VCFs for CBMR

##How to install
clone this repository with:
`git clone git@github.com:vincentrose88/SNPextractor.git`

which will create the folder `SNPextractor` from where you will run the program

##How to use

Run the extractor.sh script with the following flags and arguments (also shown if --help or -h is typed):
 * `--snp` File containing the list of SNPs that needs to extracted with the following requirement:
  * Each line has variant name, chromosome and position, seperated by `\t` (tab) (see example file). 
  * **NB: Must be on same build as imputed data (commonly b37, but do check)**. 
 * `--vcf` Path to vcf containing the HRC imputed data. 
  * For example: `/emc/cbmr/data/imputed/decode-Nov09-sanger-hwe10e5/`
 * `--type` Specify which type of data for each SNP is needed. **Only one is allowed**. 
  * Options are: DOSAGE, GENOTYPE or LIKELIHOOD. Default: DOSAGE
 * Either of the following but not both.
  * `--cohort` Prefix for a specific cohort. Extract all individuals which matches this prefix. For example `57x` matches the CIMT cohort
  * `--ind` File containing a list of individuals you want extracted. One individuals particid on each line, using "\n".
 * Optional parameters
  * `--ld` LD cut-off (for R^2) for excluding variants which are in LD with eachothers in the dataset. Only the variant with the highest imputation quality (info score) will be kept.
  * `--out` output filename. The suffix `.csv` will be added regardless. Default is: `SNPsExtracted.csv`

## Example
Extract dosage information for SNPs not in LD (R^2 < 0.5) for all inter99 individuals from Decode dataset:

`./extractor --snp SNPlist.example --vcf /emc/cbmr/data/imputed/decode-Nov09-sanger-hwe10e5/ --type DOSAGE --cohort 1000 --ld 0.5 --out inter99.decode.extract`

