#=============================================================================#
# BMI - bmi5 results
# 09/11/2019
# 
# MALES
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())

# annotations
fh_annotations <- fread("~/data/Annotations/crc_gwas_indep_signals_140_EasyStrata_LDBased.tsv")

# Rsq estimate based on alt allele probability
# filtered maf > 0.001, Rsq >= 0.8
rsq_filter <- readRDS("~/data/Rsq_Estimate/FIGI_RsqEstimate_chrALL.rds")

# results
gxe_all <- do.call(rbind, lapply(list.files(path = "~/data/Results/bmi5_male/", full.names = T, pattern = "FIGI_GxESet_bmi5_age_pc10_studygxe"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = ":"),
         chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE) %>% 
  filter(!duplicated(ID)) # small issue with gxescan package

gxe <- gxe_all %>% 
  filter(ID %in% rsq_filter$id)

# gxe_noGWAS <- filter(gxe, !SNP %in% fh_annotations$SNP)

# output chiSqGxE results for LD clumping
# calculate_pval <- function(data, statistic, df) {
#   data$P <- pchisq(data[,statistic], df = df, lower.tail = F)
#   data
# }
# 
# for(chr in 1:22) {
#   out <- calculate_pval(gxe, 'chiSqGxE', df = 1) %>%
#     filter(Chromosome == chr) %>% mutate(SNP = ID) %>%
#     dplyr::select(SNP, P)
#   write.table(out, file = paste0("/media/work/tmp/Plink_bmi5_male_ldclump_chiSqGxE_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
# }





# LD Clump Results

## chiSqGxE
tmp <- do.call(rbind, lapply(list.files("~/data/Results/bmi5_male/clump/", full.names = T, pattern = "clumped"), fread, stringsAsFactors = F))
gxe_clump_chiSqGxE <- gxe %>% 
  filter(ID %in% tmp$SNP)
rm(tmp)


#-----------------------------------------------------------------------------#
# QQ and Manhattan Plots ----
#-----------------------------------------------------------------------------#
plot_exposure <- "bmi5"
plot_covariates <- c("age_ref_imp", "study_gxe", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

# Marginal G Results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_MALE")
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_MALE")

# GxE results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGxE', df = 1, filename_suffix = "_MALE")
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGxE', df = 1, filename_suffix = "_MALE")

# 2DF results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq2df', df = 2, filename_suffix = "_MALE")
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq2df', df = 2, filename_suffix = "_MALE")

# 3DF results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq3df', df = 3, filename_suffix = "_MALE")
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq3df', df = 3, filename_suffix = "_MALE")

# GE, Case, Control ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGE', df = 1, filename_suffix = "_MALE")
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGE', df = 1, filename_suffix = "_MALE")
create_qqplot_ge(gxe, plot_exposure, plot_covariates, stat = 'chiSqControl', df = 1, filename_suffix = "_MALE")
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqControl', df = 1, filename_suffix = "_MALE")
create_qqplot_ge(gxe, plot_exposure, plot_covariates, stat = 'chiSqCase', df = 1, filename_suffix = "_MALE")
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqCase', df = 1, filename_suffix = "_MALE")



# D|G 2-step Kooperberg ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_MALE")

# G|E 2-step Murcray ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_MALE")

# EDGE 2-step Gauderman ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_MALE")


#-----------------------------------------------------------------------------#
# Two-Step after LD Clump ----
# chiSqGxE
#-----------------------------------------------------------------------------#

# D|G 2-step Kooperberg ----
gxe_twostep <- format_twostep_data(dat = gxe_clump_chiSqGxE, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = '_clump_chiSqGxE_MALE')

# G|E 2-step Murcray ----
gxe_twostep <- format_twostep_data(dat = gxe_clump_chiSqGxE, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_clump_chiSqGxE_MALE")

# EDGE 2-step Gauderman ----
gxe_twostep <- format_twostep_data(dat = gxe_clump_chiSqGxE, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_clump_chiSqGxE_MALE")
