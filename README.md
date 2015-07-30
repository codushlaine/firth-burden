# firth-burden
For binary phenotypes with one or more covariates, run a gene-based burden test

Summary
-------
0. Read in a count of variants, e.g. loss-of-function variants, for one or more genes
0. Read in a phenotype file with one or more phenotypes
0. Read in a covariate file
0. Read in a number of the minimum number of individuals with the phenotype required before running the test
0. Read in a number of the minimum number of cased carriers required before running the test
 
Examples
-------
example.gene_burden - LoF burden (counts of LoFs) per gene for 3 genes
example.phenotypes - 5 phenotypes, header is an ICD9 code in this case  
example.covariates - appropriate covariates one would use, here: SEX AGE AGESQ PC1 PC2 PC3 PC4 

Usage
-----
```Rscript run_firth_burden.R example.gene_burden example.phenotypes example.covariates 20 0 ```

Output
-----
"run_firth_burden.log" - simple log file: phenotype tested, phenotype index, time
"results/" - directory with a subdirectory for each phenotype. Within each subdirectory, each gene result is given.

Example of output for a gene:

```X8.6_pheno 109 30376 0 11 0 0.00037 2.42826 -2.43135 4.48818 11.33917 0.087918 88.95951 0.2215 2 2 DPM1_gene ```

Fields are:

PHENOTYPE / N_CASES / N_CONTROLS / CASE_CARRIERS / CONTROL_CARRIERS / CASE_FREQ / CONTROL_FREQ / BETA / BETA_CIL / BETA_CIU / OR / OR_CIL / OR_CIU / P / N_LOF / N_IND_LOF / GENE_MASK


Requirements in R
-----
library(reshape)
library(stringi)
library(logistf)

