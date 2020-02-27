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
filename <- args[2] # ex: FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan
covariates <- c(args[3:length(args)])
annotation_file <- 'gwas_140_ld_annotation.txt'

# set working directory
# important to do because figifs functions outputs into folders based on current position (e.g. figures)
setwd(exposure)

gxe <- readRDS(paste0("~/data/results/", exposure, "/processed/", filename, "_results.rds"))

#-----------------------------------------------------------------------------#
# Two-step methods
#-----------------------------------------------------------------------------#

# incorporate GWAS set main effects statistics to the two-step methods
gwas_results <- readRDS("~/data/results/gwas/processed/FIGI_v2.3_gwasset_basic_covars_gxescan_results_short.rds")

gxe_twostep_tmp <- gxe %>%
  dplyr::select(-betaG, -chiSqG, -chiSqEDGE, -chiSq3df) %>%
  inner_join(gwas_results, 'ID') %>%
  mutate(chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE)

# D|G 2-step Kooperberg
gxe_twostep <- format_twostep_data(dat = gxe_twostep_tmp, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')
# create_twostep_weighted_plot(gxe_twostep, exposure = 'asp_ref', covars = 'test', sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')

# G|E 2-step Murcray
gxe_twostep <- format_twostep_data(dat = gxe_twostep_tmp, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE')

# EDGE 2-step Gauderman
gxe_twostep <- format_twostep_data(dat = gxe_twostep_tmp, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE')

#-----------------------------------------------------------------------------#
# Two-step methods - expectation based
#-----------------------------------------------------------------------------#

# D|G 2-step Kooperberg ----
gxe_twostep_expectation <- format_twostep_data_expectation(dat = gxe_twostep_tmp, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot_expectation(gxe_twostep_expectation, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_expectation_based")

# G|E 2-step Murcray ----
gxe_twostep_expectation <- format_twostep_data_expectation(dat = gxe_twostep_tmp, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot_expectation(gxe_twostep_expectation, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_expectation_based")

# EDGE 2-step Gauderman ----
gxe_twostep_expectation <- format_twostep_data_expectation(dat = gxe_twostep_tmp, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot_expectation(gxe_twostep_expectation, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_expectation_based")


#-----------------------------------------------------------------------------#
# Two-step methods - ld clump
#-----------------------------------------------------------------------------#

tmp <- do.call(rbind, lapply(list.files(paste0("~/data/results/", exposure, "/clump/"), full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))
# tmp <- do.call(rbind, lapply(list.files(paste0("~/data/results/", 'asp_ref', "/clump/"), full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))
gxe_chiSqGxE_ld <- gxe %>%
  filter(ID %in% tmp$SNP)
rm(tmp)

# incorporate GWAS set main effects statistics to the two-step methods
gwas_results <- readRDS("~/data/results/gwas/processed/FIGI_v2.3_gwasset_basic_covars_gxescan_results_short.rds")

gxe_twostep_tmp <- gxe_chiSqGxE_ld %>%
  dplyr::select(-betaG, -chiSqG, -chiSqEDGE, -chiSq3df) %>%
  inner_join(gwas_results, 'ID') %>%
  mutate(chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE)


# D|G 2-step Kooperberg
gxe_twostep <- format_twostep_data(dat = gxe_twostep_tmp, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_chiSqGxE_ld_clump")

# G|E 2-step Murcray
gxe_twostep <- format_twostep_data(dat = gxe_twostep_tmp, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_chiSqGxE_ld_clump")

# EDGE 2-step Gauderman
gxe_twostep <- format_twostep_data(dat = gxe_twostep_tmp, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_chiSqGxE_ld_clump")

