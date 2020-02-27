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

gxetmp <- readRDS(paste0("~/data/Results/", exposure, "/processed/", filename, "_results.rds"))
gwashits <- readRDS("~/data/Annotations/gwas_140_snp_ld_and_region_based_to_exclude_from_gxe.rds")

#-----------------------------------------------------------------------------#
# remove gwas hits + SNPs r2 > 0.2
# keep hits that are 'novel' within same GWAS regions but not in LD with
# reported top hit
#-----------------------------------------------------------------------------#
annotation <- fread("~/data/Annotations/gwas_140_ld_annotation.txt") %>% 
  dplyr::mutate(SNP = paste0(Chr, ":", Pos))

gxe_nogwas <- gxetmp %>% 
  dplyr::filter(!SNP %in% annotation$SNP)

create_manhattanplot(gxe_nogwas, 'gwas', c('age_ref_imp', 'sex', 'study_gxe', 'PC1', 'PC2', 'PC3'), stat = 'chiSqG', annotation_file = 'temp_annotation_ver2.txt', df = 1, filename_suffix = "_NO_GWAS_HIT")




#-----------------------------------------------------------------------------#
# Two-step methods
#-----------------------------------------------------------------------------#
gxe <- gxetmp %>% 
  dplyr::filter(!SNP %in% gwashits$SNP)

# D|G 2-step Kooperberg ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_NO_GWAS_REGION")

# G|E 2-step Murcray ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_NO_GWAS_REGION")

# EDGE 2-step Gauderman ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_NO_GWAS_REGION")




#-----------------------------------------------------------------------------#
# Two-step methods LD CLUMPED
#-----------------------------------------------------------------------------#
tmp <- do.call(rbind, lapply(list.files(paste0("~/data/Results/", exposure, "/clump_controls/"), full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))
gxe_chiSqGxE_ld <- gxe %>% 
  filter(ID %in% tmp$SNP)
rm(tmp)


# D|G 2-step Kooperberg ----
gxe_twostep <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_LD_Clumped_NO_GWAS_REGION")

# G|E 2-step Murcray ----
gxe_twostep <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_LD_Clumped_NO_GWAS_REGION")

# EDGE 2-step Gauderman ----
gxe_twostep <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = exposure, covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_LD_Clumped_NO_GWAS_REGION")



