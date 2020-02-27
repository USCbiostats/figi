#====================================================
# FIGI GxEScan run on Merged Files
#
# John test - GetVCFInfo
#====================================================
library(BinaryDosage)

args <- commandArgs(trailingOnly=T)
chr <- args[1]

# write BDose using *.vcf.gz files
vcfInfo <- GetVCFInfo(paste0("/auto/pmd-02/figi/HRC_VCF_SampleRename/UKB1_chr", chr, ".vcf.gz"), gzipped = TRUE, index = FALSE)
bdfilenames <- c(paste0("/staging/dvc/andreeki/BD/UKB1_chr", chr, ".bdose"))
bdoptions <- character(0)

x <- BinaryDosage:::WriteBinaryDosageHeader(4, 2, bdfilenames, vcfInfo, bdoptions)
bdInfo <- BinaryDosage:::ReadBinaryDosageHeader(bdfilenames)
bdWriteInfo <- BinaryDosage:::AllocateBinaryDosageWriteMemory(bdInfo)
BinaryDosage:::VCFApply(vcfInfo, BinaryDosage:::WriteBinaryDosageData, bdWriteInfo)
BinaryDosage:::WriteBinaryDosageIndices(bdWriteInfo)
