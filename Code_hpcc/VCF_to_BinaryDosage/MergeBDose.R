#====================================================
# Testing the updated merge function
#
#====================================================
library(BinaryDosage)

args <- commandArgs(trailingOnly=T)
chr <- args[1]

etal <- c(paste0("/auto/pmd-02/figi/HRC_BDose/FIGI_noUKB_chr", chr, ".bdose"),
          paste0("/staging/dvc/andreeki/BD/ukbiobank_chr", chr, ".bdose"))

mergedFile <- paste0("/auto/pmd-02/figi/HRC_BDose/FIGI_chr", chr, ".bdose")

MergeBD(mergedFile, etal)
