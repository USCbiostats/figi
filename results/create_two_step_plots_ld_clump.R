#=============================================================================#
# FIGI
#
# two-step plots - ld clumped results
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())

# exposure = 'aspirin'
# covariates = c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3')
# filename <- 'FIGI_v2.3_gxeset_aspirin_basic_covars_gxescan'

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
# Two-step methods - ld clump
# 
# based on chiSqGxE ld clumping
#-----------------------------------------------------------------------------#
tmp <- do.call(rbind, lapply(list.files(paste0("~/data/results/", exposure, "/clump/"), full.names = T, pattern = "*chiSqGxE_ldclump.clumped"), fread, stringsAsFactors = F))

gxe_chiSqGxE_ld <- gxe %>%
  filter(SNP %in% tmp$SNP)
rm(tmp)

# D|G 2-step Kooperberg
gxe_twostep <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_chiSqGxE_ld_clump")
out <- filter(gxe_twostep, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqG_", exposure, "_", paste0(covariates, collapse = "_"), "_chiSqGxE_ld_clump.rds"), version = 2)

# G|E 2-step Murcray
gxe_twostep <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_chiSqGxE_ld_clump")
out <- filter(gxe_twostep, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqGE_", exposure, "_", paste0(covariates, collapse = "_"), "_chiSqGxE_ld_clump.rds"), version = 2)

# EDGE 2-step Gauderman
gxe_twostep <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_chiSqGxE_ld_clump")
out <- filter(gxe_twostep, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqEDGE_", exposure, "_", paste0(covariates, collapse = "_"), "_chiSqGxE_ld_clump.rds"), version = 2)



#-----------------------------------------------------------------------------#
# Two-step methods - ld clump
# 
# based on chiSqG ld clumping
#-----------------------------------------------------------------------------#
tmp <- do.call(rbind, lapply(list.files(paste0("~/data/results/", exposure, "/clump/"), full.names = T, pattern = "*chiSqG_ldclump.clumped"), fread, stringsAsFactors = F))

gxe_chiSqG_ld <- gxe %>%
  filter(SNP %in% tmp$SNP)
rm(tmp)

# D|G 2-step Kooperberg
gxe_twostep <- format_twostep_data(dat = gxe_chiSqG_ld, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_chiSqG_ld_clump")
out <- filter(gxe_twostep, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqG_", exposure, "_", paste0(covariates, collapse = "_"), "_chiSqG_ld_clump.rds"), version = 2)

# G|E 2-step Murcray
gxe_twostep <- format_twostep_data(dat = gxe_chiSqG_ld, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_chiSqG_ld_clump")
out <- filter(gxe_twostep, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqGE_", exposure, "_", paste0(covariates, collapse = "_"), "_chiSqG_ld_clump.rds"), version = 2)

# EDGE 2-step Gauderman
gxe_twostep <- format_twostep_data(dat = gxe_chiSqG_ld, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_chiSqG_ld_clump")
out <- filter(gxe_twostep, step2p < wt)
saveRDS(out, file = paste0("files/twostep_wht_chiSqEDGE_", exposure, "_", paste0(covariates, collapse = "_"), "_chiSqG_ld_clump.rds"), version = 2)




# # incorporate GWAS set main effects statistics to the two-step methods
# gwas_results <- readRDS("~/data/results/gwas/processed/FIGI_v2.3_gwasset_basic_covars_gxescan_results_short.rds")
# 
# gxe_twostep_tmp <- gxe_chiSqGxE_ld %>%
#   dplyr::select(-betaG, -chiSqG, -chiSqEDGE, -chiSq3df) %>%
#   inner_join(gwas_results, 'ID') %>%
#   mutate(chiSqEDGE = chiSqG + chiSqGE,
#          chiSq3df = chiSqG + chiSqGxE + chiSqGE)

