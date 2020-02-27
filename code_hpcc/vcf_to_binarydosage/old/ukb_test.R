#====================================================
# FIGI GxEScan run on Merged Files
#
# Script to actually run GxEScanR
# (called from the job script generated)
#====================================================
library(GxEScanR)

# args
args <- commandArgs(trailingOnly=T)
chr <- args[1]

# components
#cov <- readRDS('FIGI_GxESet_GxEScanR_N86639_20181105.rds')
#cov <- readRDS('FIGI_GxESet_sex_87914.rds')
cov <- readRDS('FIGI_GxESet_UKB_test_sex_14885.rds')

#bdose_info <- GxEScanR::GetBinaryDosageInfo(paste0('/staging/dvc/andreeki/bdose/ukbiobank_chr', chr, '.bdose'))
bdose_info <- GxEScanR::GetBinaryDosageInfo(paste0('/staging/dvc/andreeki/bdose/ukbiobank_chr', chr, '.bdose'))

# run GxEScan
start_time <- Sys.time()
GxEScan(cov, bdose_info, paste0('ukb_original_results_sex_14885_chr', chr, ".out"))
end_time <- Sys.time()

message(paste("GxEScanR runtime: ", difftime(end_time, start_time, units = 'mins')))