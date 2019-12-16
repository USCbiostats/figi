#=============================================================================#
# FIGI - create plots
#
# create two-step weighted hypothesis testing plot for LD Clumped results
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())

# exposure = 'asp_ref'
# covariates = c('age_ref_imp', 'sex', 'study_gxe', 'PC1', 'PC2', 'PC3')
# filename <- 'FIGI_GxESet_asp_ref_age_sex_pc3_studygxe_72820'

args <- commandArgs(trailingOnly=T)
exposure <- args[1] # ex: asp_ref
filename <- args[2] # ex: FIGI_GxESet_asp_ref_age_sex_pc3_studygxe_72820
covariates <- c(args[3:length(args)])
annotation_file <- 'temp_annotation_ver2.txt'

# set working directory
# important to do because figifs functions outputs into folders based on current position (e.g. figures)
setwd(exposure)

gxe <- readRDS(paste0("~/data/Results/", exposure, "/processed/", filename, "_results.rds"))

tmp <- do.call(rbind, lapply(list.files(paste0("~/data/Results/", exposure, "/clump_controls/"), full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))
gxe_chiSqGxE_ld <- gxe %>% 
  filter(ID %in% tmp$SNP)
rm(tmp)


#-----------------------------------------------------------------------------#
# Two-step methods
#-----------------------------------------------------------------------------#

# incorporate GWAS set main effects statistics to the two-step methods
gwas_results <- readRDS("~/data/Results/gwas/processed/FIGI_GwasSet_age_sex_pc3_studygxe_124702_results_short.rds") %>% 
  dplyr::select(-betaGE, -chiSqGE)

gxe_twostep_tmp <- gxe_chiSqGxE_ld %>% 
  dplyr::select(-betaG, -chiSqG, -chiSqEDGE, -chiSq3df) %>% 
  inner_join(gwas_results, 'ID') %>% 
  mutate(chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE)


# D|G 2-step Kooperberg ----
gxe_twostep <- format_twostep_data(dat = gxe_twostep_tmp, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_LD_Clumped")

# G|E 2-step Murcray ----
gxe_twostep <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_LD_Clumped")

# EDGE 2-step Gauderman ----
gxe_twostep <- format_twostep_data(dat = gxe_twostep_tmp, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_LD_Clumped")

