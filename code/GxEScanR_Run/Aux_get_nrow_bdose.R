# helper script to get number of samples/SNPs in a binary dosage files

library(GxEScanR) 

setwd("/auto/pmd-02/figi/HRC_BDose")

for(chr in c(1:22) {
    bdose <- GxEScanR::GetBinaryDosageInfo(paste0("axiom_mecc_cfr_ky_chr", chr))
    nrow(bdose[[7]])
}