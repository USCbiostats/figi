#====================================================
# Testing the updated merge function
#
#====================================================
library(BinaryDosage)

#bd_file1 <- "/staging/dvc/andreeki/bdose/ukbiobank_chr5.bdose"
#bd_file1 <- "/auto/pmd-02/figi/HRC_BDose/ccfr_omni_chr5.bdose"
#bd_file2 <- "/auto/pmd-02/figi/HRC_BDose/mecc_chr5.bdose"

args <- commandArgs(trailingOnly=T)
chr <- args[1]

etal <- c(paste0("/staging/dvc/andreeki/bdose/FIGI_CORECT_chr", chr, ".bdose"),
          paste0("/staging/dvc/andreeki/bdose/FIGI_GECCO_chr", chr, ".bdose"))



#ukb <- "/staging/dvc/andreeki/bdose/ukbiobank_chr5.bdose"

#ukb_etal <- c(etal, ukb)

mergedFile <- paste0("/auto/pmd-02/figi/HRC_BDose/FIGI_chr", chr, ".bdose")

MergeBD(mergedFile, etal)