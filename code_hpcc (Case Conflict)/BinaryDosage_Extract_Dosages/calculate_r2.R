#=============================================================================#
# FIGI - Extract Dosages from BinaryDosage files
# 11/03/2019
#
# Take a long vector of index positions (CHROMOSOME SPECIFIC)
#   split vector into chunks, then loop (or apply) over them to extract dosages
#   output chunked dosages into rds files
#   combine rds files into a single data.frame for convenience
#=============================================================================#
library(BinaryDosage)
library(tidyverse)
library(data.table)

# always take two arguments - chr and snplist (vector of index positions)
args <- commandArgs(trailingOnly=T)
filename <- args[1] # just a string of the file name containing SNP index positions in BinaryDosage files
tophit <- as.character(args[2])


# vcfids for controls
figi_controls_vcfid <- readRDS("/staging/dvc/andreeki/gwas_annotation_ld/FIGI_controls_vcfid_73601.rds")

# combine temporary files into a single data.frame. more convenient
file_list <- list.files("/staging/dvc/andreeki/gwas_annotation_ld", pattern = filename, full.names = T)

merge.all <- function(x, y) {
    merge(x, y, all = T, by = "vcfid")
}

read.all <- function(x) {
    readRDS(x) %>%
        dplyr::filter(vcfid %in% figi_controls_vcfid$vcfid)
}

dosages <- Reduce(merge.all, lapply(file_list, read.all)) %>%
    dplyr::select(-vcfid)

# pull out top hit
# calculate r2 with every snp in data.frame

#tophit_dosage <- dosages[, tophit]
tophit_dosage <- dplyr::select(dosages, tophit)

out <- apply(dosages, 2, function(x) cor(x, tophit_dosage)^2)

saveRDS(out, file = paste0("/staging/dvc/andreeki/gwas_annotation_ld/", filename, "_r2.rds"))



