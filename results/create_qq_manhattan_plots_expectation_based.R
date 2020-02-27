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
annotation_file <- 'gwas_140_ld_annotation.txt'

# set working directory
# important to do because figifs functions outputs into folders based on current position (e.g. figures)
setwd(exposure)

gxe <- readRDS(paste0("~/data/results/", exposure, "/processed/", filename, "_results.rds"))
gxe <- readRDS("~/data/results/asp_ref/processed/FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan_results.rds")


#-----------------------------------------------------------------------------#
# Two-step methods - expectation based bin assignment
#-----------------------------------------------------------------------------#

# incorporate GWAS set main effects statistics to the two-step methods
gwas_results <- readRDS("~/data/results/gwas/processed/FIGI_v2.3_gwasset_basic_covars_gxescan_results_short.rds")

gxe_twostep_tmp <- gxe %>% 
  dplyr::select(-betaG, -chiSqG, -chiSqEDGE, -chiSq3df) %>% 
  inner_join(gwas_results, 'ID') %>% 
  mutate(chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE)


# D|G 2-step Kooperberg ----
gxe_twostep_expectation <- format_twostep_data_expectation(dat = gxe_twostep_tmp, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot_expectation(gxe_twostep_expectation, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_Expectation_Based")

# G|E 2-step Murcray ----
gxe_twostep_expectation <- format_twostep_data_expectation(dat = gxe_twostep_tmp, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot_expectation(gxe_twostep_expectation, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_Expectation_Based")

# EDGE 2-step Gauderman ----
gxe_twostep_expectation <- format_twostep_data_expectation(dat = gxe_twostep_tmp, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot_expectation(gxe_twostep_expectation, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_Expectation_Based")



# ---- delete after cleaning ---- #
results <- filter(gxe_twostep_expectation, step2p < wt) 

# of these results, which fall in regions previously annotated...
annot <- fread("~/data/Annotations/temp_annotation_ver2.txt") %>% 
  mutate(SNP = paste0(Chr, ":", Pos))

results_new <- anti_join(results, annot, 'SNP')

results_known <- anti_join(results, results_new, 'SNP') %>% 
  filter(Chromosome != 5)

## some familiar loci, and some brand new one

chr <- filter(annot, Chr == 6) %>% 
  arrange(Pos)


results <- filter(gxe_twostep_expectation, step2p < wt) %>% 
  mutate(locusgroup <- kmeans(log(MapInfo))[[1]])



