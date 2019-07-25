#=============================================================================#
# Folate - folate_dietqc2 Results
# 05/18/2019
# 
# Generate Plots
# Implement 2-step methods
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(figifs)
rm(list = ls())

# annotations
fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"))

# Rsq Filter - re-estimate using alt allele probability (see BinaryDosage forked package)
rsq_filter <- readRDS("/media/work/data/Rsq_Estimates_AA/rsq_aa_genomewide_filter_pt8.rds")

# Rsq Filter - weighted average (by imputation batch sample size)
# rsq_filter <- readRDS("~/data/HRC_InfoFile_Merged/HRC_WtAvgRsq_HRCRefPanelMAFs.rds") %>% 
#   filter(Rsq_avg > 0.8)

# Read in results
# calculate chiSqEDGE and chiSq3df statistics
gxe_all <- do.call(rbind, lapply(list.files(path = "~/data/Results/folate_totqc2/", full.names = T, pattern = "results_GxE_folate_totqc2_sex_age_pc3_energytot_studygxe_54084_binCovF_chr"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = ":"),
         chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE) %>% 
  filter(!duplicated(ID)) # small issue with gxescan package

gxe <- gxe_all %>% 
  dplyr::filter(ID %in% rsq_filter$SNP)

# output results for LD clumping
# this need to be done after filtering for Rsq values
calculate_pval <- function(data, statistic, df) {
  data$P <- pchisq(data[,statistic], df = df, lower.tail = F)
  data
}
for(chr in 1:22) {
  out <- calculate_pval(gxe, 'chiSqGxE', df = 1) %>%
    filter(Chromosome == chr) %>% mutate(SNP = ID) %>%
    dplyr::select(SNP, P)
  write.table(out, file = paste0("/media/work/data/tmp_files/Plink_ldclump_chiSqGxE_folate_totqc2_sex_age_pc3_energytot_studygxe_54084_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
}


# LD Clump Results
gxe_chiSqGxE_ldclumped <- do.call(rbind, lapply(list.files("~/data/Results/folate_totqc2/ldclump/", full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))

gxe_chiSqGxE_ld <- gxe %>% 
  filter(ID %in% gxe_chiSqGxE_ldclumped$SNP)
rm(gxe_chiSqGxE_ldclumped)



#-----------------------------------------------------------------------------#
# QQ and Manhattan Plots ----
#-----------------------------------------------------------------------------#
plot_exposure <- "folate_totqc2"
plot_covariates <- c("age_ref_imp", "sex", "study_gxe", "PC1", "PC2", "PC3", "energytot")

# Marginal G Results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1)

# GxE results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGxE', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGxE', df = 1)

# 2DF results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq2df', df = 2)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq2df', df = 2)

# 3DF results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq3df', df = 3)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq3df', df = 3)

# GE, Case, Control ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGE', df = 1)
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqControl', df = 1)
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqCase', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqCase', df = 1)



# D|G 2-step Kooperberg ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')

# G|E 2-step Murcray ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE')

# EDGE 2-step Gauderman ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE')



#-----------------------------------------------------------------------------#
# chiSqGxE ld clumped output
#-----------------------------------------------------------------------------#
# D|G 2-step Kooperberg ----
gxe_twostep_ld <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep_ld, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')

# G|E 2-step Murcray ----
gxe_twostep_ld <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep_ld, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE')

# EDGE 2-step Gauderman ----
gxe_twostep_ld <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep_ld, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE')







