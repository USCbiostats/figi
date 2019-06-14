# FIGI PCA  

Calculate PCs using combined FIGI data  

## Random Sample of 30,000 Imputed SNPs  

- SNP_Sample_30K.R  
  - MAF > 0.05, Rsq > 0.99 (based on *corect_oncoarray* imputation batch)  
  - Filter SNPs overlapping between HRC and 1KGP  

## Summary  

- j01.sh - extract markers, convert to plink format  
- j02.sh - merge chromosomes  
- j03.sh - merge imputation batches  
  - file name: **ALL_Merged_PCA**
- j04.sh - extract markers (KGP files)
- j05.sh - merge FIGI + KGP  
  - file name: **ALL_Merged_PCA_KGP**  
- j06.sh - run PCA (plink2 --approx)  