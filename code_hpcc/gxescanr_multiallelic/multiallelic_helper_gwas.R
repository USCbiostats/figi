library(tidyverse)
library(data.table)
library(BinaryDosage)
library(GxEScanR)

wrap <- function(chr, exposure) {
    snps <- readRDS(paste0("/staging/dvc/andreeki/gxescanr_multiallelic/hrc_multiallelic_sites_index_chr", chr, ".rds"))
    geno <- readRDS(paste0("/auto/pmd-02/figi/HRC_BDose/FIGI_snpid_fix_chr", chr, ".rds"))
    pheno <- readRDS(paste0("/staging/dvc/andreeki/gxescanr_multiallelic/", exposure, "/FIGI_v2.3_gxeset_", exposure, "_basic_covars_gxescan.rds"))
    GxEScan(subjectData=pheno, geneticData=geno, binCov = F, outFile=paste0("/staging/dvc/andreeki/gxescanr_multiallelic/", exposure, "/FIGI_v2.3_gxeset_", exposure, "_basic_covars_gxescan_multiallelic_chr", chr, ".out"), skipFile=paste0("/staging/dvc/andreeki/gxescanr_multiallelic/", exposure, "/FIGI_v2.3_gxeset_", exposure, "_basic_covars_gxescan_multiallelic_skipped_chr", chr, ".out"), popminMaf = 0, sampleminMaf=0.01, snps = snps)
}

