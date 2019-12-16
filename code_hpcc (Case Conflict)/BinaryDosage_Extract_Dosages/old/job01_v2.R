#========================================================#
# obtain allele dosage from FIGI BDose files
# using index position
#========================================================#
library(BinaryDosage)

args <- commandArgs(trailingOnly=T)
chr <- args[1]
snplist <- args[2]

# run function
snpsToRead <- readRDS(paste0(snplist, ".rds"))
bdose <- readRDS(paste0("/auto/pmd-02/figi/HRC_BDose/FIGI_chr", chr, ".rds")) # FIGI_chr.rds index file
snpsbd <- GetSNPValues(bdose, snpsToRead, geneProb = F)

out <- data.frame(snpsbd)
out$vcfid <- rownames(out)

saveRDS(out, file = paste0(snplist, "_out.rds"))
