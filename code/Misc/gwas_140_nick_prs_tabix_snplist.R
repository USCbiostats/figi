#=============================================================================#
# gecco gwas 140 (jeroen natgen) 
# output effects + risk allele for nick
# use meta-analysis results Flora sent you (125k)
# this script outputs list of SNPs to use with tabix - more convenient)
#
#
# use this with tabix indexed results
# "MarginalMeta_HRC_EUR_only_Results.tsv.gz"
#=============================================================================# 
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())


# Jeroen GWAS hits
# important - EXCLUDE 6:105966894 (does not appear to be in HRC reference panel)
gwas <-fread("~/Dropbox/FIGI/Documentation/Annotations_GWASHits_20190913/crc_gwas_indep_signals_050719_temp.tsv")[,-c(16,17)] %>% mutate(SNP = paste(CHROM, POSITION, sep = ":")) %>% 
  arrange(SNP)

any(duplicated(gwas$RSID))


out <- gwas %>% 
  dplyr::mutate(tabix = paste0(CHROM, ":", POSITION-1, "-", POSITION)) %>% 
  dplyr::mutate(chrom = CHROM,
                # chrom = paste0("chr", CHROM),
                chromStart = POSITION-1,
                chromEnd = POSITION)

write.table(out[, c("chrom", "chromStart", "chromEnd")], file = "~/gecco_125k_gwas_eur/gwas_140_nick_prs.bed", quote = F, row.names = F, col.names = F, sep = '\t')
