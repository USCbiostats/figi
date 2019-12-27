#=============================================================================#
# FIGI - create plots
#
# 1) QQ + Manhattan Plots for GxEScanR outputs
# one major change is to make all figure names consistent so you can 
# generalize the rmarkdown report generation. I'm tired. 
#
# 2) Create plots for LD Clumped results in a separate script
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
# setwd(paste0("~/Dropbox/FIGI/results/", exposure))

# output from preprocess_gxescanr_results.R
gxe <- readRDS(paste0("~/data/results/", exposure, "/processed/", filename, "_results.rds"))

#-----------------------------------------------------------------------------#
# QQ and Manhattan Plots
#-----------------------------------------------------------------------------#

# Marginal G Results ----
create_qqplot(gxe, exposure, covariates, stat = 'chiSqG', df = 1)
create_manhattanplot(gxe, exposure, covariates, stat = 'chiSqG', annotation_file = annotation_file, df = 1)

# GxE results ----
create_qqplot(gxe, exposure, covariates, stat = 'chiSqGxE', df = 1)
create_manhattanplot(gxe, exposure, covariates, stat = 'chiSqGxE', annotation_file = annotation_file,df = 1)

# 2DF results ----
create_qqplot(gxe, exposure, covariates, stat = 'chiSq2df', df = 2)
create_manhattanplot(gxe, exposure, covariates, stat = 'chiSq2df', annotation_file = annotation_file,df = 2)

# 3DF results ----
create_qqplot(gxe, exposure, covariates, stat = 'chiSq3df', df = 3)
create_manhattanplot(gxe, exposure, covariates, stat = 'chiSq3df', annotation_file = annotation_file,df = 3)

# GE ----
create_qqplot(gxe, exposure, covariates, stat = 'chiSqGE', df = 1)
create_manhattanplot(gxe, exposure, covariates, stat = 'chiSqGE', annotation_file = annotation_file, df = 1)

# GE Among Controls ----
create_qqplot_ge(gxe, exposure, covariates, stat = 'chiSqControl', df = 1)
create_manhattanplot(gxe, exposure, covariates, stat = 'chiSqControl', annotation_file = annotation_file, df = 1)

# GE Case Only
create_qqplot_ge(gxe, exposure, covariates, stat = 'chiSqCase', df = 1)
create_manhattanplot(gxe, exposure, covariates, stat = 'chiSqCase', annotation_file = annotation_file, df = 1)


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

# D|G 2-step Kooperberg ----
gxe_twostep <- format_twostep_data(dat = gxe_twostep_tmp, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')

# G|E 2-step Murcray ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE')

# EDGE 2-step Gauderman ----
gxe_twostep <- format_twostep_data(dat = gxe_twostep_tmp, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE')
