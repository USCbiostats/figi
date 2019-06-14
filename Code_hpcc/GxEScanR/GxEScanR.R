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
GxEOnly(figiCov, figiGene, paste0("/staging/dvc/andreeki/figi_gxescan_gxeonly/results_GxE_asp_ref_sex_age_pc10_studygxe_noUKB_58109_chr", chr, ".out"), 0.01, 0.05)
Sys.time()







class(bdose) <- genetic-file-info