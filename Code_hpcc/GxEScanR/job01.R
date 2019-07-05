#=============================================================================#
# GxEScanR
# FIGI GxESet (N = 102792)
# 
# exposure: asp_ref (aspirin/nsaids at ref)
# model: outcome ~ age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe + G*asp_ref
# complete case N = 72820
#=============================================================================#
library(GxEScanR)
args <- commandArgs(trailingOnly=T)
chr <- args[1]


wdir <- "/auto/pmd-02/figi/HRC_BDose"
covfile <- "/staging/dvc/andreeki/GxE/FIGI_GxESet_folate_dietqc2_sex_age_pc3_energytot_studygxe_52447.rds"
bdosefile <- paste0("FIGI_chr", chr, ".rds") 

outFile <-      paste0("/staging/dvc/andreeki/GxE/results_GxE_folate_dietqc2_sex_age_pc3_energytot_studygxe_52447_binCovF_chr", chr, ".out")
outFile_skip <- paste0("/staging/dvc/andreeki/GxE/results_GxE_folate_dietqc2_sex_age_pc3_energytot_studygxe_52447_binCovF_SKIPPED_chr", chr, ".out")



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
