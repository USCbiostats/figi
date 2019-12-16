#=============================================================================#
# Calculate Alternate Allele Frequencies
# function in BinaryDosage package
#=============================================================================#
library(BinaryDosage)

setwd("/auto/pmd-02/figi/HRC_BDose")
args <- commandArgs(trailingOnly=T)
chr <- args[1]

# Sample list (vector of vcfids)
samples <- readRDS("/staging/dvc/andreeki/AAF/FIGI_GWASSet_vcfid_N_138015.rds")

# Read BinaryDosage information file (index)
figiGene <- readRDS(paste0("FIGI_chr", chr, ".rds"))

out <- GetAlternateAlleleFrequencies(figiGene, 1L:nrow(figiGene[[11]]), samples)
saveRDS(out, file = paste0("/staging/dvc/andreeki/AAF/FIGI_AAF_chr", chr, ".rds"))
