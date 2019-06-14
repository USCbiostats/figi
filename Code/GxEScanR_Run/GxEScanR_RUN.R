#====================================================
# FIGI GxEScan run on Merged Files
#
# Script to actually run GxEScanR 
# (called from the job script generated)
#====================================================
library(tidyverse)
library(data.table)
library(GxEScanR)

# ARGS
args <- commandArgs(trailingOnly=T)
filename <- args[1]
chr <- args[2]
covar <- args[3]
chunk <- args[4]

# locations
jobid <- paste0(Sys.Date(), "_", covar, "/") # should always be the same honestly
jobid <- paste0(covar, "/") # should always be the same honestly
work_dir <- "/staging/dvc/andreeki/figi/"
jobs_dir <- paste0("/staging/dvc/andreeki/figi/jobs/", jobid)
logs_dir <- paste0("/staging/dvc/andreeki/figi/logs/", jobid)
results_dir <- paste0("/staging/dvc/andreeki/figi/results/", jobid)
bdose_dir <- "/auto/pmd-02/figi/HRC_BDose/"
setwd(work_dir)

# components
bdose_info <- GxEScanR::GetBinaryDosageInfo(paste0(bdose_dir, filename, "_chr", chr, ".bdose"))
cov <- read.table(paste0("~/", covar, ".txt"), stringsAsFactors = F, header = T)
snps <- readRDS(paste0(jobs_dir, filename, "_chr", chr, "_chunk_", chunk, ".rds"))

# run GxEScan
start_time <- Sys.time()
GxEScan(cov, bdose_info, paste0(results_dir, filename, "_chr", chr,"_chunk_", chunk, "_", covar, ".out"), snps = snps)
end_time <- Sys.time()

start_time
end_time
end_time - start_time
