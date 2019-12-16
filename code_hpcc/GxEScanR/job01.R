#=============================================================================#
# GxEScanR
# FIGI GxESet (N = 102792)
#=============================================================================#
library(GxEScanR)
args <- commandArgs(trailingOnly=T)
chr <- args[1]
covf <- args[2]

wdir <- "/auto/pmd-02/figi/HRC_BDose"
covfile <- paste0(covf, ".rds")
bdosefile <- paste0("FIGI_chr", chr, ".rds") 

outFile <-      paste0(covf, "_chr", chr, ".out")
outFile_skip <- paste0(covf, "_SKIPPED_chr", chr, ".out")


#-----------------------------------------------------------------------------#
# set directory
setwd(wdir)

# Read in the covariate info
figiCov <- readRDS(covfile)

# Read in the information of the binary dosage file
figiGene <- readRDS(bdosefile)
class(figiGene) <- "genetic-file-info" # (temporary, john will fix)


#numSNPs <- nrow(figiGene[[11]])



# Fit the models
Sys.time()
#GxEScan(figiCov, figiGene, outFile=outFile, skipFile=outFile_skip, popminMaf=0.01, sampleminMaf=0.01, binCov=F, snps = 2678934:numSNPs)
GxEScan(figiCov, figiGene, outFile=outFile, skipFile=outFile_skip, popminMaf=0.01, sampleminMaf=0.01, binCov=F)
Sys.time()
