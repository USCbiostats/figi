#=============================================================================#
# Calculate Rsq in pooled imputed data (based on minimac3 formula)
# function in BinaryDosage package
# 
# Probably best to use the GWAS set (~ 138K individuals)... right? 
#=============================================================================#
library(BinaryDosage)

setwd("/auto/pmd-02/figi/HRC_BDose")
args <- commandArgs(trailingOnly=T)
chr <- args[1]

# vcfid vector (GWAS or GxE sets)
samples <- readRDS("/staging/dvc/andreeki/Rsq/FIGI_GWASSet_vcfid_N_138015.rds")

# binarydosage info file (indexed)
figiGene <- readRDS(paste0("FIGI_chr", chr, ".rds"))

out <- GetRsqEstimateMinimac(figiGene, 1L:nrow(figiGene[[11]]), samples)
saveRDS(out, file = paste0("/staging/dvc/andreeki/Rsq/FIGI_RsqEstimate_chr", chr, ".rds"))
