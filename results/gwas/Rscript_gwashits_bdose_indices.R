#=============================================================================#
# 09/10/2019
#
# Extract GWAS top hits (N = 140)
# for exclusion + PRS analysis
#
#=============================================================================#
library(tidyverse)
library(data.table)
rm(list = ls())

# single file, no LD based annotation, just top hits
annot <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_061918.tsv")[, -c('LL_95PCT_CI_BETA', 'UL_95PCT_CI_BETA')] %>% 
  mutate(ID = gsub("/", ":", gsub("_", ":", VARIANT))) %>% 
  arrange(CHROM, POSITION)


# how about a function that takes snps in the format chr:bp:ref:alt, and then just outputs a vector (after looping/applying) of index positions in a figi_chrNN.rds object. so all you have to provide is a vector of SNPs to extract without having to do this again and again on my end

# get_snp_index <- function(snp) {
#   chr <- strsplit(snp, "[:]")[[1]][1]
#   snpinfo <- readRDS(paste0("~/data/BinaryDosage_InfoFile/FIGI_chr", chr, ".rds"))[[11]]
#   snpinfo$ID <- paste(snpinfo$SNPID, snpinfo$Reference, snpinfo$Alternate, sep = ":")
#   as.integer(match(snp, snpinfo$ID))
# }
# 
# get_snp_index("20:33213196:A:C")


# on second thought, no matter what, you'll have to send files to HPC, might as well give it vectors of indices right away.. 
# so the way I was doing it was fine

annot_chroms <- as.integer(names(table(annot$CHROM)))

for(chr in annot_chroms) {
  figi <- readRDS(paste0("~/data/BinaryDosage_InfoFile/FIGI_chr", chr, ".rds"))[[11]] %>% 
    mutate(ID = paste(Chromosome, Location, Reference, Alternate, sep = ":"))
  saveRDS(which(figi$ID %in% annot$ID), file = paste0("files/GWAS_hits_indices_chr", chr, ".rds"), version = 2)
}

x <- readRDS("files/GWAS_hits_indices_chr19.rds")

chr <- 1
figi <- readRDS(paste0("~/data/BinaryDosage_InfoFile/FIGI_chr", chr, ".rds"))[[11]] %>% 
  mutate(ID = paste(Chromosome, Location, Reference, Alternate, sep = ":"))
wtf <- which(figi$ID %in% annot$ID)


wtf <- readRDS("~/Dropbox/FIGI/Results/gwas/files/GWAS_hits_indices_chr1.rds")
