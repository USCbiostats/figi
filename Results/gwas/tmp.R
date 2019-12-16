#=============================================================================#
# FIGI - create plots
# after filtering out gwas hit regions (I created a file...)
#
# by filtering gwas hits (main effects), you're really just trying to gauge how
# different the two-step method plots will perform
#
# also filter by GE statistic.. eventually
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())

# set working directory
# important to do because figifs functions outputs into folders based on current position (e.g. figures)
setwd("~/Dropbox/FIGI/Results/gwas")

gxetmp <- readRDS("~/data/Results/gwas/processed/FIGI_GwasSet_age_sex_pc3_studygxe_124702_results.rds")
gwashits <- readRDS("~/data/Annotations/gwas_140_snp_ld_and_region_based_to_exclude_from_gxe.rds")

gxe <- gxetmp %>% 
  dplyr::filter(!SNP %in% gwashits$SNP)


# Marginal G Results ----
create_manhattanplot(gxe, "gwas", c("age", "sex", "study_gxe", "PC1", "PC2", "PC3"), stat = 'chiSqG', annotation_file = 'temp_annotation_ver2.txt', df = 1, filename_suffix = "_NO_GWAS_HITS")



