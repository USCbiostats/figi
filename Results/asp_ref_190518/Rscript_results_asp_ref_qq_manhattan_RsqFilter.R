#=============================================================================#
# NSAIDS - asp_ref Results
# 05/18/2019
# 
# Generate Plots
# Implement 2-step methods
#
# let's see what things look like after filtering by Rsq
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

# Rsq Filter (vector of IDs to keep)
rsq_filter <- readRDS("~/data/HRC_InfoFile_Merged/HRC_WtAvgRsq_HRCRefPanelMAFs.rds") %>% 
  filter(Rsq_avg > 0.8)

# results
gxe_all <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref_190518/", full.names = T, pattern = "results_GxE_asp_ref_sex_age_pc3_studygxe_72820_binCovT_chr"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = ":"),
         chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE) %>% 
  filter(!duplicated(ID))

gxe <- gxe_all %>% 
  filter(ID %in% rsq_filter$ID)


# you want to make sure NOT to drop annotations from manhattan plots. let's make sure there's good overlap...
# 90 OF THE MARKERS ARE IN MY RESULTS, MISSING 5, 2 of them were previously reported
# annotation_original <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_061918.tsv")[, -c('LL_95PCT_CI_BETA', 'UL_95PCT_CI_BETA')] %>% 
#   filter(SUGGESTIVE=="no") %>% 
#   separate(VARIANT, into = c("SNP", "alleles"), sep = "_")
# 
# check <- inner_join(gxe, annotation_original, by = "SNP")
# check <- anti_join(annotation_original, gxe, by = "SNP")



# output chiSqGxE results for LD clumping
# gxe <- gxe %>% 
#   dplyr::mutate(P = pchisq(chiSqGxE, df = 1, lower.tail = F))
# 
# for(chr in 1:22) {
#   out <- gxe %>% 
#     filter(Chromosome == chr) %>% mutate(SNP = ID) %>%
#     dplyr::select(SNP, P)
#   write.table(out, file = paste0("/media/work/tmp/LDclump_chiSqGxE_asp_ref_age_ref_imp_sex_study_gxe_PC1-3_N_72820_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
# }


#-----------------------------------------------------------------------------#
# QQ and Manhattan Plots ----
#-----------------------------------------------------------------------------#
plot_exposure <- "asp_ref"
plot_covariates <- c("age_ref_imp", "sex", "study_gxe", "PC1", "PC2", "PC3")

# Marginal G Results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1)

# GxE results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGxE', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGxE', df = 1)

# 2DF results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq2df', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq2df', df = 1)

# 3DF results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq3df', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq3df', df = 1)

# GE, Case, Control ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGE', df = 1)
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqControl', df = 1)
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqCase', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqCase', df = 1)



# D|G 2-step Kooperberg ----
gxe_twostep <- format_2step_data(data = gxe, 'chiSqG', 5, 0.05)
create_2step_weighted_plot(gxe_twostep, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')

# G|E 2-step Murcray ----
gxe_twostep <- format_2step_data(data = gxe, 'chiSqGE', 5, 0.05)
create_2step_weighted_plot(gxe_twostep, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE')

# EDGE 2-step Gauderman ----
gxe_twostep <- format_2step_data(data = gxe, 'chiSqEDGE', 5, 0.05)
create_2step_weighted_plot(gxe_twostep, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE')





