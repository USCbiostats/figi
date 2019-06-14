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
cov <- readRDS('ukbiobank_chr5.bdoseInfo.rds')
bdose_info <- readRDS('ukbiobank_chr5.bdoseInfo.rds')

# run GxEScan
start_time <- Sys.time()
GxEScan(cov, bdose_info, paste0('results_N86639_20181105_chr', chr, ".out"))
end_time <- Sys.time()

message(paste("GxEScanR runtime: ", difftime(end_time, start_time, units = 'mins')))