# using gxeonly version gxescanr (01/17/2019)
library(GxEScanR)
setwd("/auto/pmd-02/figi/HRC_BDose/FIGI_NeedFixAsians")
args <- commandArgs(trailingOnly=T)
chr <- args[1]

# Read in the covariate info
figiCov <- readRDS("/staging/dvc/andreeki/figi_gxescan_gxeonly/FIGI_GxESet_asp_ref_sex_age_pc10_studygxe_noUKB_58109.rds")

# Read in the information of the binary dosage file
figiGene <- readRDS(paste0("FIGI_chr", chr, ".rds"))

# Fit the models
Sys.time()
GEOnly(figiCov, figiGene, paste0("/staging/dvc/andreeki/figi_gxescan_gxeonly/results_GEOnly_asp_ref_sex_age_pc10_studygxe_noUKB_58109_chr", chr, ".out"), 0.01, 0.05)
Sys.time()