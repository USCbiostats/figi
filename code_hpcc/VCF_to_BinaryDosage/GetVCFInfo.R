#====================================================
# FIGI GxEScan run on Merged Files
#
# John test - GetVCFInfo
#====================================================
library(BinaryDosage)

# args
args <- commandArgs(trailingOnly=T)
vcfFile <- args[1]

# components
#cov <- readRDS('FIGI_GxESet_GxEScanR_N86639_20181105.rds')
#cov <- readRDS('FIGI_GxESet_sex_87914.rds')
#bdose_info <- GxEScanR::GetBinaryDosageInfo(paste0('/auto/pmd-02/figi/HRC_BDose/FIGI_ALL_chr', chr, '.bdose'))

# run GxEScan
start_time <- Sys.time()
#vcfinfo <- GetVCFInfo(vcfFile=vcfFile, reserve=3000000)
vcfinfo <- GetVCFInfo(vcfFile=vcfFile)
end_time <- Sys.time()

saveRDS(vcfinfo, file = paste0(vcfFile, ".rds"))

message(paste("GetVCFInfo runtime: ", difftime(end_time, start_time, units = 'mins')))