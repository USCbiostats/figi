#=============================================================================#
# GxEScanR
# FIGI GxESet (N = 102792)
#=============================================================================#
library(BinaryDosage)
args <- commandArgs(trailingOnly=T)
chr <- args[1]

wdir <- "/auto/pmd-02/figi/HRC_BDose"
bdosefile <- paste0("FIGI_chr", chr, ".rds") 

#-----------------------------------------------------------------------------#
# set directory
setwd(wdir)

# Read in the information of the binary dosage file
figiGene <- readRDS(bdosefile)

# Fit the models
Sys.time()
out <- GetAlternateAlleleFrequencies(figiGene, SNPs = 1L:nrow(figiGene[[11]]))
saveRDS(out, file = "/staging/dvc/andreeki/test.rds")
Sys.time()
