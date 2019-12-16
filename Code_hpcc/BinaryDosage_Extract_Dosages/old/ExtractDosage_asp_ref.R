#========================================================#
# obtain allele dosage from FIGI BDose files
# using index position
#========================================================#
library(BinaryDosage)
args <- commandArgs(trailingOnly=T)
chr <- args[1]

setwd("/auto/pmd-02/figi/HRC_BDose/FIGI_NeedFixAsians/")

snpsToRead <- readRDS(paste0("/staging/dvc/andreeki/gxescanr/GxEScanR_asp_ref_top1262_bdose_SNP_chr", chr, ".rds"))
bd <- readRDS(paste0("FIGI_chr", chr, ".rds"))

snpsbd <- GetSNPValues(bd, snpsToRead, geneProb = F)
saveRDS(snpsbd, file = paste0("/staging/dvc/andreeki/gxescanr/GxEScanR_asp_ref_extract_chr", chr, ".rds"))