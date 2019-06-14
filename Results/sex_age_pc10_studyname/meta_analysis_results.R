#========================================================#
# Manhattan and QQ Plots
# (qqman)
#========================================================#
library(tidyverse)
library(data.table)
library(qqman)

rm(list = ls())

# jeroen meta-analysis results
joint_chr15 <- fread("~/results/jeroen_meta/colorectal_cancer_gwas_summary_statistics/METAANALYSIS_JOINT_CHR15_1_SORTED.tbl") %>% 
  arrange(`P-value`)



# top hits from stephanie
annot <- fread("~/bin/EasyStrata/crc_gwas_125k_indep_signals_061918.tsv")

# 
# [, -c('LL_95PCT_CI_BETA', 'UL_95PCT_CI_BETA')] %>%
#   filter(SUGGESTIVE == "no") %>%
#   dplyr::rename(Chr = CHROM,
#                 Pos = POSITION) %>%
#   mutate(Colour = "gold2")
# write.table(annot, file = "~/bin/EasyStrata/GWAS_TopHits_95.txt", quote = F, row.names = F, sep = '\t')

# original results
results <- do.call(rbind, lapply(list.files(full.names = T), fread, stringsAsFactors = F))
names(results)


