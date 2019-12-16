#=============================================================================#
# GxEScanR
# FIGI GxESet (N = 102792)
# 
# exposure: bmi5
# model: outcome ~ age_ref_imp + sex + PC1-10 + study_gxe + bmi5 + G + bmi5*G
# complete case N = 88,324
#=============================================================================#
library(GxEScanR)
args <- commandArgs(trailingOnly=T)
chr <- args[1]


wdir <- "/auto/pmd-02/figi/HRC_BDose"
covfile <- "/staging/dvc/andreeki/GxE/bmi5/FIGI_GxESet_bmi5_sex_age_pc10_studygxe_88324.rds"
bdosefile <- paste0("FIGI_chr", chr, ".rds") 

outFile <-      paste0("/staging/dvc/andreeki/GxE/bmi5/FIGI_GxESet_bmi5_sex_age_pc10_studygxe_88324_chr", chr, ".out")
outFile_skip <- paste0("/staging/dvc/andreeki/GxE/bmi5/FIGI_GxESet_bmi5_sex_age_pc10_studygxe_88324_SKIPPED_chr", chr, ".out")



#-----------------------------------------------------------------------------#
# set directory
setwd(wdir)

# Read in the covariate info
figiCov <- readRDS(covfile)

# Read in the information of the binary dosage file
figiGene <- readRDS(bdosefile)
class(figiGene) <- "genetic-file-info" # (temporary, john will fix)

# Fit the models
Sys.time()
GxEScan(figiCov, figiGene, outFile=outFile, skipFile=outFile_skip, popminMaf=0.01, sampleminMaf=0.01, binCov=F)
Sys.time()
