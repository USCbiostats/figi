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
gxe_all <- do.call(rbind, lapply(list.files(path = "~/data/Results/folate_dietqc2/", full.names = T, pattern = "results_GxE_folate_dietqc2_sex_age_pc3_energytot_studygxe_52447_binCovF_chr"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = ":"),
         chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE) %>% 
  filter(!duplicated(ID)) # small issue with gxescan package

gxe <- gxe_all %>% 
  dplyr::filter(ID %in% rsq_filter$SNP)

# output results for LD clumping
# this need to be done after filtering for Rsq values
# calculate_pval <- function(data, statistic, df) {
#   data$P <- pchisq(data[,statistic], df = df, lower.tail = F)
#   data
# }
# for(chr in 1:22) {
#   out <- calculate_pval(gxe, 'chiSqGxE', df = 1) %>%
#     filter(Chromosome == chr) %>% mutate(SNP = ID) %>%
#     dplyr::select(SNP, P)
#   write.table(out, file = paste0("/media/work/data/tmp_files/Plink_clump_chiSqGxE_folate_dietqc2_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
# }


# LD Clump Results
gxe_chiSqGxE_ldclumped <- do.call(rbind, lapply(list.files("~/data/Results/folate_dietqc2/ld_clump/", full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))

gxe_chiSqGxE_ld <- gxe %>% 
  filter(ID %in% gxe_chiSqGxE_ldclumped$SNP)
rm(gxe_chiSqGxE_ldclumped)




#-----------------------------------------------------------------------------#
# QQ and Manhattan Plots ----
#-----------------------------------------------------------------------------#
plot_exposure <- "folate_dietqc2"
plot_covariates <- c("age_ref_imp", "sex", "study_gxe", "PC1", "PC2", "PC3", "energytot")

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
gxe_twostep <- format_data_twostep_data(data = gxe, 'chiSqG', 5, 0.05)
create_2step_weighted_plot(gxe_twostep, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')

# G|E 2-step Murcray ----
gxe_twostep <- format_2step_data(data = gxe, 'chiSqGE', 5, 0.05)
create_2step_weighted_plot(gxe_twostep, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE')

# EDGE 2-step Gauderman ----
gxe_twostep <- format_2step_data(data = gxe, 'chiSqEDGE', 5, 0.05)
create_2step_weighted_plot(gxe_twostep, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE')






















#-----------------------------------------------------------------------------#
# Marginal G Results ----
#-----------------------------------------------------------------------------#
# QQ Plot
create_qqplot(gxe, 'chiSqG', df = 1)
# Manhattan Plot
create_manhattanplot(gxe, 'chiSqG', df = 1)

#-----------------------------------------------------------------------------#
# GxE results ----
#-----------------------------------------------------------------------------#
# QQ Plot
create_qqplot(gxe, 'chiSqGxE', df = 1)
# Manhattan Plot
create_manhattanplot(gxe, 'chiSqGxE', df = 1)

# single hit on chr6, almost significant actual peak looking on chr1
tmp <- gxe %>% 
  filter(chiSqGxE > 20)

#-----------------------------------------------------------------------------#
# 2DF results ----
#-----------------------------------------------------------------------------#
# QQ Plot
create_qqplot(gxe, 'chiSq2df', df = 2)
# Manhattan Plot
create_manhattanplot(gxe, 'chiSq2df', df = 2)


#-----------------------------------------------------------------------------#
# 3DF results ----
#-----------------------------------------------------------------------------#
# QQ Plot
create_qqplot(gxe, 'chiSq3df', df = 3)
# Manhattan Plot
create_manhattanplot(gxe, 'chiSq3df', df = 3)



#-----------------------------------------------------------------------------#
# Miami Plot G vs 2DF ----
#-----------------------------------------------------------------------------#
# Miami Plot
create_miamiplot(gxe, 'chiSq2df', 'chiSqG', 2, 1)

#-----------------------------------------------------------------------------#
# GE, Case, Control ----
#-----------------------------------------------------------------------------#

# GE QQ Plot
create_qqplot(gxe, 'chiSqGE', df = 1)

# Control-Only
create_qqplot(gxe, 'chiSqControl', df = 1)

# Case-Only QQ + manhattan
create_qqplot(gxe, 'chiSqCase', df = 1)
create_manhattanplot(gxe, 'chiSqCase', df = 1)


#-----------------------------------------------------------------------------#
# D|G 2-step Kooperberg ----
#-----------------------------------------------------------------------------#
gxe_normal <- format_2step_data(data = gxe, 'chiSqG', 5, 0.05)
create_2step_weighted_plot(gxe_normal, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')

# ld clumped results
gxe_twostep_ld <- format_2step_data(data = gxe_ld, 'chiSqG', 5, 0.05)
create_2step_weighted_plot(gxe_twostep_ld, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')

#-----------------------------------------------------------------------------#
# G|E 2-step Murcray ----
#-----------------------------------------------------------------------------#
gxe_normal <- format_2step_data(data = gxe, 'chiSqGE', 5, 0.05)
create_2step_weighted_plot(gxe_normal, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE')

# ld clumped results
gxe_twostep_ld <- format_2step_data(data = gxe_ld, 'chiSqGE', 5, 0.05)
create_2step_weighted_plot(gxe_twostep_ld, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE')


#-----------------------------------------------------------------------------#
# EDGE 2-step Gauderman ----
#-----------------------------------------------------------------------------#
gxe_normal <- format_2step_data(data = gxe, 'chiSqEDGE', 5, 0.05)
create_2step_weighted_plot(gxe_normal, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE')

# ld clumped results
gxe_twostep_ld <- format_2step_data(data = gxe_ld, 'chiSqEDGE', 5, 0.05)
create_2step_weighted_plot(gxe_twostep_ld, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE')


#-----------------------------------------------------------------------------#
# 3DF MIAMI (vs G) ----
#-----------------------------------------------------------------------------#
# Miami Plot
create_miamiplot(gxe, 'chiSq3df', 'chiSqG', 3, 1)


#-----------------------------------------------------------------------------#
# output for Prioritypruner ----
#-----------------------------------------------------------------------------#

# Required Columns:
#   name - Name of the SNP (e.g., rs222).
# chr - Name of the chromosome (e.g., 1, chr1). Chromosome X must be denoted as either 'X', 'chrX' or '23'. Chromosome Y and MT are not supported.
# pos - Physical position of the SNP.
# a1 - First allele of the SNP.
# a2 - Second allele of the SNP.
# p - P-value or other prioritization metric between 0 and 1. This is used for prioritizing the selection of SNPs, where lower numbers are prioritized.
# forceSelect - Flag indicating if the SNP should be selected (kept) regardless of its LD with other selected SNPs or other filtering criteria specified, such as MAF or design score (1=true, 0=false).
# designScore - Design score of the SNP (any positive real number). Can be filled in with a constant value (e.g., 1) if unknown.


# let's start with p values based on D|G statistic.. 

# didn't filter by UKB... don't filter here

results_G <- gxe %>%
  rename(name = SNP, 
         chr = Chromosome, 
         pos = Location,
         a1 = Reference,
         a2 = Alternate) %>% 
  mutate(p = pchisq(chiSqG, df = 1, lower.tail = F),
         forceSelect = 0,
         designScore = 1) %>%
  # filter(!ID %in% ukb_filter$ID,
  #        chr == 22) %>% 
  filter(chr == 22) %>% 
  dplyr::select(name, chr, pos, a1, a2, p, forceSelect, designScore)

write.table(results_G, file = "~/test.txt", quote = F, row.names = F)

# UGH need to update sex information on tfam file
cov <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc10_studygxe_72145_GLM.rds")

fam <- fread("~/final.tfam") %>% 
  inner_join(cov)

# also keep in mind multiallelics...
# simply removed them all. 
# code got lost because of rstudio crash, but easily write again. 



+-
# results

wtfff <- fread("~/wtfff.results")


