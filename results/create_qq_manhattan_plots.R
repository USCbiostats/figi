#=============================================================================#
# FIGI - Create QQ and Manhattan Plots
# 
# QQ + Manhattan Plots for GxEScanR outputs
#
# I COULD generate 2df/3df using gwas marginal but i don't think it's 
# necessary. Leave it for two-step plots
#
#
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())

# exposure = 'asp_ref'
# covariates = c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3')
# filename <- 'FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan'

args <- commandArgs(trailingOnly=T)
exposure <- args[1] # ex: asp_ref
filename <- args[2] # ex: FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan
covariates <- c(args[3:length(args)])
annotation_file <- 'gwas_140_ld_annotation_new.txt'

# set working directory
# important to do because figifs functions outputs into folders based on current position (e.g. figures)
setwd(paste0("~/Dropbox/FIGI/Results/", exposure))

# output from preprocess_gxescanr_results.R
gxe <- readRDS(paste0("~/data/results/", exposure, "/processed/", filename, "_results.rds"))
# gxe <- readRDS(paste0("~/data/results/asp_ref/processed/FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan_results.rds"))


#-----------------------------------------------------------------------------#
# QQ and Manhattan Plots ----
#-----------------------------------------------------------------------------#

# Marginal G Results
create_qqplot(gxe, exposure, covariates, stat = 'chiSqG', df = 1)
create_manhattanplot(gxe, exposure, covariates, stat = 'chiSqG', annotation_file = annotation_file, df = 1)

# GxE results
create_qqplot(gxe, exposure, covariates, stat = 'chiSqGxE', df = 1)
create_manhattanplot(gxe, exposure, covariates, stat = 'chiSqGxE', annotation_file = annotation_file,df = 1)

# 2DF results
create_qqplot(gxe, exposure, covariates, stat = 'chiSq2df', df = 2)
create_manhattanplot(gxe, exposure, covariates, stat = 'chiSq2df', annotation_file = annotation_file,df = 2)

# 3DF results
create_qqplot(gxe, exposure, covariates, stat = 'chiSq3df', df = 3)
create_manhattanplot(gxe, exposure, covariates, stat = 'chiSq3df', annotation_file = annotation_file,df = 3)

# GE
create_qqplot(gxe, exposure, covariates, stat = 'chiSqGE', df = 1)
create_manhattanplot(gxe, exposure, covariates, stat = 'chiSqGE', annotation_file = annotation_file, df = 1)

# GE Among Controls
create_qqplot_ge(gxe, exposure, covariates, stat = 'chiSqControl', df = 1)
create_manhattanplot(gxe, exposure, covariates, stat = 'chiSqControl', annotation_file = annotation_file, df = 1)

# GE Case Only
create_qqplot_ge(gxe, exposure, covariates, stat = 'chiSqCase', df = 1)
create_manhattanplot(gxe, exposure, covariates, stat = 'chiSqCase', annotation_file = annotation_file, df = 1)


#-----------------------------------------------------------------------------#
# QQ and Manhattan Plots, no GWAS hits ----
#-----------------------------------------------------------------------------#

snps_to_exclude <- readRDS("~/data/Annotations/gwas_140_snp_ld_and_region_based_to_exclude_from_gxe.rds")

gxe_nogwas <- gxe %>%
  mutate(chr_bp = paste0(Chromosome, ":", Location)) %>%
  dplyr::filter(!chr_bp %in% snps_to_exclude$SNP)


# Marginal G Results
# create_qqplot(gxe, exposure, covariates, stat = 'chiSqG', df = 1)
create_manhattanplot(gxe_nogwas, exposure, covariates, stat = 'chiSqG', annotation_file = annotation_file, df = 1, filename_suffix = "_no_gwas")

# 2DF results ----
# create_qqplot(gxe_nogwas, exposure, covariates, stat = 'chiSq2df', df = 2, filename_suffix = "_no_gwas")
create_manhattanplot(gxe_nogwas, exposure, covariates, stat = 'chiSq2df', annotation_file = annotation_file,df = 2, filename_suffix = "_no_gwas")

# 3DF results ----
# create_qqplot(gxe_nogwas, exposure, covariates, stat = 'chiSq3df', df = 3, filename_suffix = "_no_gwas")
create_manhattanplot(gxe_nogwas, exposure, covariates, stat = 'chiSq3df', annotation_file = annotation_file,df = 3, filename_suffix = "_no_gwas")
