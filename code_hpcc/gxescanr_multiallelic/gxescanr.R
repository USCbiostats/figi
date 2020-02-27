#=============================================================================#
# GxEScanR
# 
# ------ Notes: ------ #
# - MAF < 1%
#=============================================================================#
library(GxEScanR)
args <- commandArgs(trailingOnly=T)
chr <- args[1]
cov <- args[2]

bdose_directory <- "/auto/pmd-02/figi/HRC_BDose"
covariate_file <- paste0(cov, ".rds")
bdose_file <- paste0("FIGI_snpid_fix_chr", chr, ".rds") 

outFile <-      paste0(cov, "_chr", chr, ".out")
outFile_skip <- paste0(cov, "_skipped_chr", chr, ".out")


#-----------------------------------------------------------------------------#
# set directory
setwd(bdose_directory)

# Read in the covariate info
figiCov <- readRDS(covariate_file)

# Read in the information of the binary dosage file
figiGene <- readRDS(bdose_file)
class(figiGene) <- "genetic-file-info" # (temporary, john will fix)

# Fit the models
Sys.time()
#GxEScan(figiCov, figiGene, outFile=outFile, skipFile=outFile_skip, popminMaf=0.01, sampleminMaf=0.01, binCov=F, snps = 2678934:numSNPs)
#GxEScan(figiCov, figiGene, outFile=outFile, skipFile=outFile_skip, popminMaf=0.01, sampleminMaf=0.01, binCov=F)
GxEScan(figiCov, figiGene, outFile=outFile, skipFile=outFile_skip, popminMaf=0, sampleminMaf=0.01, binCov=F)
Sys.time()
