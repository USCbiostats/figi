#=============================================================================#
# NSAIDS - aspirin Results
# 06/01/2019
# 
# Filter by Weighted Average Rsq (until another solution is devised)
# Generate Plots
# Implement 2-step methods
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
rm(list = ls())
source("~/Dropbox/FIGI/Code/Functions/GxEScan_PostHoc_Analyses.R")


# annotations
# (right now, I only use this df to ensure that these markers stay in plot after filtering by P < 0.05)
fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"))

# Rsq Filter (vector of IDs to keep) --- will need to recreate this at some point
rsq_filter <- readRDS("~/data/HRC_InfoFile_Merged/HRC_WtAvgRsq_HRCRefPanelMAFs.rds") %>% 
  filter(Rsq_avg > 0.8)

# results
gxe_all <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/aspirin_190528/", full.names = T, pattern = "results_GxE_aspirin_sex_age_pc3_studygxe_66485_binCovT_chr"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = ":"),
         chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE) %>% 
  filter(!duplicated(ID)) # small issue with gxescan package

gxe <- gxe_all %>% 
  filter(ID %in% rsq_filter$ID)


# output chiSqGxE results for LD clumping
calculate_pval <- function(data, statistic, df) {
  data$P <- pchisq(data[,statistic], df = df, lower.tail = F)
  data
}

for(chr in 1:22) {
  out <- calculate_pval(gxe, 'chiSqGxE', df = 1) %>%
    filter(Chromosome == chr) %>% mutate(SNP = ID) %>%
    dplyr::select(SNP, P)
  write.table(out, file = paste0("/media/work/tmp/Plink_aspirin_ldclump_chiSqGxE_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')

}



# LD Clump Results
gxe_ldclumped <- do.call(rbind, lapply(list.files("~/data/Results/NSAIDS/aspirin_190528/aspirin_chiSqGxE_ldclump/", full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))

gxe_ld <- gxe %>% 
  filter(ID %in% gxe_ldclumped$SNP)


# ----- testing ------ #
# output results for LD clumping chiSqG (to compare manhattan plots)
for(chr in 1:22) {
  out <- calculate_pval(gxe, 'chiSqG', df = 1) %>%
    filter(Chromosome == chr) %>% mutate(SNP = ID) %>%
    dplyr::select(SNP, P)
  write.table(out, file = paste0("/media/work/tmp/Plink_aspirin_ldclump_chiSqG_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
}

gxe_chiSqG_ldclumped <- do.call(rbind, lapply(list.files("~/data/Results/NSAIDS/aspirin_190528/aspirin_chiSqG_ldclump/", full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))

gxe_ld <- gxe %>% 
  filter(ID %in% gxe_chiSqG_ldclumped$SNP)



#-------------------------------------#
# convenience (labels, titles etc)
# VERY IMPORTANT TO DEFINE THESE !!!
# (CLEAN ENVIRONMENT BEFORE RUNNING!)
#-------------------------------------#
global_N <- unique(gxe$Subjects)
global_covs <- c("age_ref_imp", "sex", "study_gxe", "PC1-3")
global_E <- "aspirin"


#-----------------------------------------------------------------------------#
# Marginal G Results ----
#-----------------------------------------------------------------------------#
# QQ Plot
create_qqplot(gxe, 'chiSqG', df = 1)
# Manhattan Plot
create_manhattanplot(gxe, 'chiSqG', df = 1)

# Manhattan Plot
create_manhattanplot(gxe_ld, 'chiSqG', df = 1)



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

# ------ Clumped results ------ #
x <- gxe %>% 
  filter(ID %in% gxe_ldclumped$SNP)

create_manhattanplot(x, 'chiSqGxE', df = 1)


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
tmp <- format_2step_data(data = gxe, 'chiSqG', 5, 0.05)
create_2step_weighted_plot(tmp, sizeBin0 = 5, alpha = 0.05, binsToPlot = 15, statistic = 'chiSqG')


# test plotting only ld clumped results
tmp <- format_2step_data(data = gxe_ld, 'chiSqG', 5, 0.05)
create_2step_weighted_plot(tmp, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')

#-----------------------------------------------------------------------------#
# G|E 2-step Murcray ----
#-----------------------------------------------------------------------------#
tmp <- format_2step_data(data = gxe, 'chiSqGE', 5, 0.05)
create_2step_weighted_plot(tmp, sizeBin0 = 5, alpha = 0.05, binsToPlot = 15, statistic = 'chiSqGE')




# test plotting only ld clumped results
tmp <- format_2step_data(data = gxe_ld, 'chiSqGE', 5, 0.05)
create_2step_weighted_plot(tmp, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE')




#-----------------------------------------------------------------------------#
# EDGE 2-step Gauderman ----
#-----------------------------------------------------------------------------#
tmp <- format_2step_data(data = gxe, 'chiSqEDGE', 5, 0.05)
create_2step_weighted_plot(tmp, sizeBin0 = 5, alpha = 0.05, binsToPlot = 15, statistic = 'chiSqEDGE')


# test plotting only ld clumped results
tmp <- format_2step_data(data = gxe_ld, 'chiSqEDGE', 5, 0.05)
create_2step_weighted_plot(tmp, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE')


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

  # results
  
  wtfff <- fread("~/wtfff.results")



