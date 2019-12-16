#=========================================================
# FIGI GxEScan run on Merged Files
# 
# Create slurm job scripts and rds objects for submission
#=========================================================
library(tidyverse)
library(data.table)
library(GxEScanR)

args <- commandArgs(trailingOnly=T)
filename <- args[1]
chr <- args[2]
covar <- args[3] 
# testing
#filename <- "FIGI_ALL"; chr <- 22; covar <- "FIGI_Covariates_09272018_age_pc_platform_sex"

# locations
jobid <- paste0(Sys.Date(), "_", covar, "/")
jobid <- paste0(covar, "/")

work_dir <- "/staging/dvc/andreeki/figi/"
jobs_dir <- paste0("/staging/dvc/andreeki/figi/jobs/", jobid)
logs_dir <- paste0("/staging/dvc/andreeki/figi/logs/", jobid)
results_dir <- paste0("/staging/dvc/andreeki/figi/results/", jobid)
bdose_dir <- "/auto/pmd-02/figi/HRC_BDose/"
info_dir <- "/auto/pmd-02/figi/HRC_VCF/corect_oncoarray/" # using oncoarray as standard info file
setwd(work_dir)
chunk_size <- 50000
maf_filter <- 0.05

    # (is this bad)
system(paste0("mkdir -p ", jobs_dir))
system(paste0("mkdir -p ", logs_dir))
system(paste0("mkdir -p ", results_dir))

# Remove Typed_Only markers
# Filter by MAF
bdose <- GxEScanR::GetBinaryDosageInfo(paste0(bdose_dir, filename, "_chr", chr, ".bdose"))
info <- fread(paste0("zcat ", info_dir, "chr", chr, ".info.gz")) %>%
  filter(Genotyped != 'Typed_Only')
filtr <- bind_cols(bdose[[9]], bdose[[10]]) %>% 
  mutate(MAF = 0.5 - abs(AAF - 0.5))
snps <- which(filtr$SNP %in% info$SNP & 
              filtr$MAF > maf_filter)

# save rds objects
snps_chunk <- split(snps, ceiling(seq_along(snps)/chunk_size))
lapply(seq_along(snps_chunk), function(i) saveRDS(snps_chunk[[i]], paste0(jobs_dir, filename, "_chr", chr, "_chunk_", names(snps_chunk)[[i]], ".rds")))

# create slurm job scripts
create_sbatch <- function(x) {
    jobname <- paste0(jobs_dir, filename, "_chr", chr, "_chunk_", x, ".job")
    sink(jobname)
    cat("#!/bin/bash\n")
    cat("#SBATCH --time=100:00:00\n")
    cat("#SBATCH --ntasks=1\n")
    cat("#SBATCH --mem=4GB\n")
    cat("#SBATCH --account=lc_dvc\n")
    cat("#SBATCH --partition=conti\n")
    #cat("#SBATCH --mail-type=END\n")
    cat("#SBATCH --job-name=", filename,"_chr", chr,"_chunk_", x, "\n", sep="")
    cat("#SBATCH --output=", logs_dir, filename, "_chr", chr, "_chunk_", x, ".log\n", sep="")
    cat("Rscript /staging/dvc/andreeki/figi/GxEScanR_RUN.R", filename, chr, covar, x, "\n", sep=" ")
    sink()

    # Submit to run on cluster
    #system(paste("sbatch", jobname))
}
lapply(names(snps_chunk), create_sbatch)