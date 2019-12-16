#=============================================================================#
# GxEScanR
# FIGI GxESet (N = 102792)
# 
# exposure: folate_totqc2 (total folate intake at ref)
# model: outcome ~ age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe + energytot + G*folate_totqc2
# complete case N = 54084
#=============================================================================#
library(GxEScanR)
args <- commandArgs(trailingOnly=T)
chr <- args[1]


wdir <- "/auto/pmd-02/figi/HRC_BDose"
covfile <- "/staging/dvc/andreeki/GxE/folate_totqc2/FIGI_GxESet_folate_totqc2_sex_age_pc3_energytot_studygxe_54084.rds"
bdosefile <- paste0("FIGI_chr", chr, ".rds") 

outFile <-      paste0("/staging/dvc/andreeki/GxE/folate_totqc2/results_GxE_folate_totqc2_sex_age_pc3_energytot_studygxe_54084_binCovF_chr", chr, ".out")
outFile_skip <- paste0("/staging/dvc/andreeki/GxE/folate_totqc2/results_GxE_folate_totqc2_sex_age_pc3_energytot_studygxe_54084_binCovF_SKIPPED_chr", chr, ".out")



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
