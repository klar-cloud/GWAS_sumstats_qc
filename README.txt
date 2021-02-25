Semi-automatic pipeline for preprocessing and harmonizing GWAS multiple sumstats from various sources.

Column names used in the script: BP, CHR, ID, RSID, PVAL, A1, A2, FREQ, BETA, SE, ZSCORE, NCASES, NCONTROLS, NEFF, N

Necessary software:
- R >=3.5 (modules: dplyr, data.table, optparse, rlist, ggplot2, viridis, SNPlocs.Hsapiens.dbSNP144.GRCh37)
- Python3 (modules: pandas >=1.0.0, numpy, matplotlib, argparse)

To run this pipeline you need:
- GWAS summary statistics (txt, csv) with necessary columns: BP, CHR, PVAL, A1, A2, FREQ, ZSCORE or BETA & SE
- CNAMES file - (txt, csv; space/tab-separated) file, which pairs unique names of sumstats columns and the names used in the script; Examplary file is provided with more detailed info how to create it. 	

Optional:
- REFERENCE file (txt, csv) with four columns: A1_REF, A2_REF, FREQ_REF, COORDINATE. If not provided, script will create custom reference file.


EXAMPLE QUALITY CHECK PART1
****************************

add NCASES and NCONTROLS; for continuous traits put N/2 in both variables; 
if dataset with more than one batch give a list of integers (ncases and ncontrols must have an equal length)
for datsets with efficient N column, columns with N for control&cases, 
or column with N (for continuous traits) give the approximation of values - neccessary for calculating lambda1000 in QQ plot

Functions:
- harmonize labels
- calculate efficient sample size
- create a coordinator column 
- remove ambiguous and multiallelic SNPs
- remove SNPs with low MAF
- create or update a reference panel (optional)


cd $PATH_QC/qc1
python qc1.py \
--sumstats [sumstats_path] \
--cnames [cnames_path] \
--out [sumstats_save_path] \
--ncases [ncases]\  
--ncontrols [ncontrols]\
--remove-amb \
--remove-misp \
--remove-maf [maf] \
--make-reference [ref_save_path] \
--test

EXAMPLE QUALITY CHECK PART2
****************************
Functions:
- strand-flip the mismatching alleles
- compare alleles frequency between datset and reference and removing outliers
- add RSIDs (if necessary)
- remove duplicated RSIDs
- calculate zscore (if necessary)
- infere beta and sd (if necessary)


cd $PATH_QC/qc2
Rscript qc2.R \
--sumstats [sumstats_save_path]'.qc1.csv' \
--reference [ref_save_path] \
--out [final_save_path] \
--density-plot \
--check-rsid \
--infere-beta \
--add-zscore \
--test
