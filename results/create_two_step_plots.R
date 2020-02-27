#=============================================================================#
# FIGI
#
# Create all two-step weighted hypothesis testing plots
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
filename <- args[2] # ex: FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan
covariates <- c(args[3:length(args)])
annotation_file <- 'gwas_140_ld_annotation_new.txt'

# set working directory
# important to do because figifs functions outputs into folders based on current position (e.g. figures)
setwd(paste0("~/Dropbox/FIGI/Results/", exposure))

gxe <- readRDS(paste0("~/data/results/", exposure, "/processed/", filename, "_results.rds"))


#-----------------------------------------------------------------------------#
# Two-step methods 
# 
# step1: 
# - using main effects statistics from exposure subset
# - main effects models are adjusted by exposure
#-----------------------------------------------------------------------------#


# D|G 2-step Kooperberg
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')
out <- filter(gxe_twostep, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqG_", exposure, "_", paste0(covariates, collapse = "_"), ".rds"), version = 2)

# G|E 2-step Murcray
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE')
out <- filter(gxe_twostep, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqGE_", exposure, "_", paste0(covariates, collapse = "_"), ".rds"), version = 2)

# EDGE 2-step Gauderman
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE')
out <- filter(gxe_twostep, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqEDGE_", exposure, "_", paste0(covariates, collapse = "_"), ".rds"), version = 2)





#-----------------------------------------------------------------------------#
# Two-step methods, expectation based
#
# step1: 
# - using main effects statistics from exposure subset
# - main effects models are adjusted by exposure
#-----------------------------------------------------------------------------#

# D|G 2-step Kooperberg ----
gxe_twostep_expectation <- format_twostep_data_expectation(dat = gxe, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot_expectation(gxe_twostep_expectation, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_expectation_based")
out <- filter(gxe_twostep_expectation, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqG_", exposure, "_", paste0(covariates, collapse = "_"), "_expectation_based.rds"), version = 2)

# G|E 2-step Murcray ----
gxe_twostep_expectation <- format_twostep_data_expectation(dat = gxe, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot_expectation(gxe_twostep_expectation, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_expectation_based")
out <- filter(gxe_twostep_expectation, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqGE_", exposure, "_", paste0(covariates, collapse = "_"), "_expectation_based.rds"), version = 2)

# EDGE 2-step Gauderman ----
gxe_twostep_expectation <- format_twostep_data_expectation(dat = gxe, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot_expectation(gxe_twostep_expectation, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_expectation_based")
out <- filter(gxe_twostep_expectation, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqEDGE_", exposure, "_", paste0(covariates, collapse = "_"), "_expectation_based.rds"), version = 2)





#-----------------------------------------------------------------------------#
# Two-step methods, step 1 using gwas statistics
# 
# step1: 
# - using FIGI GWAS set marginal statistics (N ~ 120K)
# - incorporate those stats into two-step methods
#-----------------------------------------------------------------------------#

# incorporate GWAS set main effects statistics to the two-step methods
gwas_results <- readRDS("~/data/results/gwas/processed/FIGI_v2.3_gwasset_basic_covars_gxescan_results_short.rds")

gxe_twostep_tmp <- gxe %>%
  dplyr::select(-betaG, -chiSqG, -chiSqEDGE, -chiSq3df) %>%
  inner_join(gwas_results, 'SNP') %>%
  mutate(chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE)


# D|G 2-step Kooperberg
gxe_twostep <- format_twostep_data(dat = gxe_twostep_tmp, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_gwas_step1")
out <- filter(gxe_twostep, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqG_", exposure, "_", paste0(covariates, collapse = "_"), "_gwas_step1.rds"), version = 2)

# G|E 2-step Murcray
gxe_twostep <- format_twostep_data(dat = gxe_twostep_tmp, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_gwas_step1")
out <- filter(gxe_twostep, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqGE_", exposure, "_", paste0(covariates, collapse = "_"), "_gwas_step1.rds"), version = 2)

# EDGE 2-step Gauderman
gxe_twostep <- format_twostep_data(dat = gxe_twostep_tmp, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_gwas_step1")
out <- filter(gxe_twostep, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqEDGE_", exposure, "_", paste0(covariates, collapse = "_"), "_gwas_step1.rds"), version = 2)



#-----------------------------------------------------------------------------#
# Two-step methods, expectation based
#
# use 
#-----------------------------------------------------------------------------#

# D|G 2-step Kooperberg ----
gxe_twostep_expectation <- format_twostep_data_expectation(dat = gxe_twostep_tmp, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot_expectation(gxe_twostep_expectation, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_gwas_step1_expectation_based")
out <- filter(gxe_twostep_expectation, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqG_", exposure, "_", paste0(covariates, collapse = "_"), "_gwas_step1_expectation_based.rds"), version = 2)

# G|E 2-step Murcray ----
gxe_twostep_expectation <- format_twostep_data_expectation(dat = gxe_twostep_tmp, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot_expectation(gxe_twostep_expectation, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_gwas_step1_expectation_based")
out <- filter(gxe_twostep_expectation, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqGE_", exposure, "_", paste0(covariates, collapse = "_"), "_gwas_step1_expectation_based.rds"), version = 2)

# EDGE 2-step Gauderman ----
gxe_twostep_expectation <- format_twostep_data_expectation(dat = gxe_twostep_tmp, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot_expectation(gxe_twostep_expectation, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_gwas_step1_expectation_based")
out <- filter(gxe_twostep_expectation, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqEDGE_", exposure, "_", paste0(covariates, collapse = "_"), "_gwas_step1_expectation_based.rds"), version = 2)


