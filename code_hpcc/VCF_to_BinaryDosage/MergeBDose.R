#====================================================
# Testing the updated merge function
#
#====================================================
library(BinaryDosage)

args <- commandArgs(trailingOnly=T)
chr <- args[1]

etal <- c(paste0("/auto/pmd-02/figi/HRC_BDose/FIGI_CORECT_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/FIGI_GECCO_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/UKB1_chr", chr, ".bdose"), 
          paste0("/auto/pmd-02/figi/HRC_BDose/UKB2_chr", chr, ".bdose"))

mergedFile <- paste0("/auto/pmd-02/figi/HRC_BDose/FIGI_chr", chr, ".bdose")

MergeBD(mergedFile, etal)
