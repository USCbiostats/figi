#=============================================================================#
# FIGI
#
# Output data.frames of GWIS significant findings
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())

args <- commandArgs(trailingOnly=T)
exposure <- args[1] # ex: asp_ref
filename <- args[2] # ex: FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan
covariates <- c(args[3:length(args)])
annotation_file <- 'gwas_140_ld_annotation.txt'

# set working directory
# important to do because figifs functions outputs into folders based on current position (e.g. figures)
setwd(paste0("~/Dropbox/FIGI/Results/", exposure))
gxe <- readRDS(paste0("~/data/results/", exposure, "/processed/", filename, "_results.rds"))


#-----------------------------------------------------------------------------#
# prepare data according to statistic, filter, and output as rds
#-----------------------------------------------------------------------------#

plot_sig_wrapper <- function(x, stat, df){
  tmp <- x %>% 
    mutate(Pval = calculate_pval(x[,stat], df)) %>% 
    dplyr::filter(Pval < 5e-8) %>%
    mutate(Pval = formatC(Pval, format = 'e', digits = 2)) %>% 
    dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval)
}

# Marginal G Results
out <- plot_sig_wrapper(gxe, 'chiSqG', df = 1)
saveRDS(out, file = paste0("files/chiSqG_", exposure, "_", paste0(covariates, collapse = "_"), ".rds"), version = 2)

# # GxE results
# out <- plot_sig_wrapper(gxe, 'chiSqGxE', df = 1)
# saveRDS(out, file = paste0("files/chiSqGxE_", exposure, "_", paste0(covariates, collapse = "_"), ".rds"), version = 2)

# 2DF results
out <- plot_sig_wrapper(gxe, 'chiSq2df', df = 2)
saveRDS(out, file = paste0("files/chiSq2df_", exposure, "_", paste0(covariates, collapse = "_"), ".rds"), version = 2)

# 3DF results
out <- plot_sig_wrapper(gxe, 'chiSq3df', df = 3)
saveRDS(out, file = paste0("files/chiSq3df_", exposure, "_", paste0(covariates, collapse = "_"), ".rds"), version = 2)

# GE
out <- plot_sig_wrapper(gxe, 'chiSqGE', df = 1)
saveRDS(out, file = paste0("files/chiSqGE_", exposure, "_", paste0(covariates, collapse = "_"), ".rds"), version = 2)

# GE Among Controls
out <- plot_sig_wrapper(gxe, 'chiSqControl', df = 1)
saveRDS(out, file = paste0("files/chiSqControl_", exposure, "_", paste0(covariates, collapse = "_"), ".rds"), version = 2)

# GE Case Only
out <- plot_sig_wrapper(gxe, 'chiSqCase', df = 1)
saveRDS(out, file = paste0("files/chiSqCase_", exposure, "_", paste0(covariates, collapse = "_"), ".rds"), version = 2)


#-----------------------------------------------------------------------------#
# prepare data according to statistic, filter, and output as rds
# for GxE - report SUGGESTIVE (5e-6)
#-----------------------------------------------------------------------------#

plot_sig_wrapper <- function(x, stat, df){
  tmp <- x %>% 
    mutate(Pval = calculate_pval(x[,stat], df)) %>% 
    dplyr::filter(Pval < 5e-6) %>%
    mutate(Pval = formatC(Pval, format = 'e', digits = 2)) %>% 
    dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval)
}

# GxE results
out <- plot_sig_wrapper(gxe, 'chiSqGxE', df = 1)
saveRDS(out, file = paste0("files/chiSqGxE_", exposure, "_", paste0(covariates, collapse = "_"), ".rds"), version = 2)

#-----------------------------------------------------------------------------#
# exclude gwas top hits 
#-----------------------------------------------------------------------------#

snps_to_exclude <- readRDS("~/data/Annotations/gwas_140_snp_ld_and_region_based_to_exclude_from_gxe.rds")

gxe_nogwas <- gxe %>%
  mutate(chr_bp = paste0(Chromosome, ":", Location)) %>%
  dplyr::filter(!chr_bp %in% snps_to_exclude$SNP)

# 2DF results
out <- plot_sig_wrapper(gxe_nogwas, 'chiSq2df', df = 2)
saveRDS(out, file = paste0("files/chiSq2df_", exposure, "_", paste0(covariates, collapse = "_"), "_no_gwas.rds"), version = 2)

# 3DF results
out <- plot_sig_wrapper(gxe_nogwas, 'chiSq3df', df = 3)
saveRDS(out, file = paste0("files/chiSq3df_", exposure, "_", paste0(covariates, collapse = "_"), "_no_gwas.rds"), version = 2)
