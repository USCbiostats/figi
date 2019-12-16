#=============================================================================#
# Calculcate r2 between top hit and SNPs within 1MB of hit
# 11/03/2019
#
# to be run after you extract dosages in chunks.. 
#
# remember youre running this separately for each top hit locus, 
# called sbatch for each of the 140 of them 
#=============================================================================#
library(BinaryDosage)
library(tidyverse)
library(data.table)

# take two arguments - chromosome and bp position of top hit
args <- commandArgs(trailingOnly=T)
chr <- args[1] # just a string of the file name containing SNP index positions in BinaryDosage files
bp <- args[2]

# vcfids for controls
figi_controls_vcfid <- readRDS("/staging/dvc/andreeki/gwas_annotation_ld/FIGI_controls_vcfid_73601.rds")


file_list <-  list.files("/staging/dvc/andreeki/gwas_annotation_ld", pattern = paste0("huygue_gwas_140_", chr, "_", bp, "_TMP"), full.names = T)
tophit <- readRDS(paste0("/staging/dvc/andreeki/gwas/huygue_gwas_140_tophitonly_", chr, "_", bp, "_TMP1.rds")) %>%
    dplyr::filter(vcfid %in% figi_controls_vcfid$vcfid)
tophit_dosage <- tophit[, 1]



# remember the rds file contains data.frame (vcfid, tophitdosage)

read_r2_wrapper <- function(x) {
    tmp <- readRDS(x) %>%
        dplyr::filter(vcfid %in% figi_controls_vcfid$vcfid) %>%
        dplyr::select(-vcfid)
    apply(tmp, 2, function(y) cor(y, tophit_dosage)^2)
}

out <- lapply(file_list, read_r2_wrapper)
out_final <- do.call(c, out)

saveRDS(out_final, file = paste0("/staging/dvc/andreeki/gwas_annotation_ld/huygue_gwas_140_", chr, "_", bp, "_r2.rds"))



