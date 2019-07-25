#=============================================================================#
# NSAIDS - aspirin results
# 06/01/2019
# 
# Filter by Weighted Average Rsq (until another solution is devised e.g. reimputation)
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
# calculate_pval <- function(data, statistic, df) {
#   data$P <- pchisq(data[,statistic], df = df, lower.tail = F)
#   data
# }
# 
# for(chr in 1:22) {
#   out <- calculate_pval(gxe, 'chiSqGxE', df = 1) %>%
#     filter(Chromosome == chr) %>% mutate(SNP = ID) %>%
#     dplyr::select(SNP, P)
#   write.table(out, file = paste0("/media/work/tmp/Plink_aspirin_ldclump_chiSqGxE_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
# 
# }



# LD Clump Results
gxe_chiSqGxE_ldclumped <- do.call(rbind, lapply(list.files("~/data/Results/NSAIDS/aspirin_190528/aspirin_chiSqGxE_ldclump/", full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))

gxe_chiSqGxE_ld <- gxe %>% 
  filter(ID %in% gxe_chiSqGxE_ldclumped$SNP)





# ----- testing ------ #
# # output results for LD clumping chiSqG (to compare manhattan plots)
# for(chr in 1:22) {
#   out <- calculate_pval(gxe, 'chiSqG', df = 1) %>%
#     filter(Chromosome == chr) %>% mutate(SNP = ID) %>%
#     dplyr::select(SNP, P)
#   write.table(out, file = paste0("/media/work/tmp/Plink_aspirin_ldclump_chiSqG_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
# }

# gxe_chiSqG_ldclumped <- do.call(rbind, lapply(list.files("~/data/Results/NSAIDS/aspirin_190528/aspirin_chiSqG_ldclump/", full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))
# 
# gxe_chiSqG_ld <- gxe %>% 
#   filter(ID %in% gxe_chiSqG_ldclumped$SNP)



#-----------------------------------------------------------------------------#
# QQ and Manhattan Plots ----
#-----------------------------------------------------------------------------#
plot_exposure <- "aspirin"
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








#-----------------------------------------------------------------------------#
# Followup interesting hits ----
#-----------------------------------------------------------------------------#
gxe_twostep_ld <- format_2step_data(data = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05) %>% 
  filter(step2p < wt)



#-----------------------------------------------------------------------------#
# locuszoom ------
# there are several regions, we should focus on one at 
# a time
#-----------------------------------------------------------------------------#

# Marginal G 

# GxE
gxe_locuszoom <- gxe %>%
  mutate(`P-value` = pchisq(chiSqGxE, df = 1, lower.tail = F),
         MarkerName = paste0("chr", SNP)) %>%
  dplyr::select(MarkerName, `P-value`)


write.table(gxe_locuszoom, file = "~/locuszoom/examples/GxEScanR_GxE_aspirin_age_ref_imp_sex_study_gxe_PC1-3_N_66485.txt", quote = F, row.names = F, sep = "\t")


#-----------------------------------------------------------------------------#
# Extract top hits dosages from binarydosage 
#-----------------------------------------------------------------------------#
# first take care of sample list (vcfid)
# might want to move that to the top of the file.. 
sample_list <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_66485_GLM.rds")[, 'vcfid']
saveRDS(sample_list, file = paste0(write_binarydosage_vcfid_filename(), ".rds"), version = 2)


# get SNPs of interest based on plots above - in this case, two-step data.table
# (don't forget in this case the sig results was from clumped file)
# do it by chromosome for now
gxe_twostep <- format_2step_data(data = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05) 
gxe_twostep_sig <- gxe_twostep[step2p < wt, ]

get_binarydosage_index(gxe_twostep_sig$ID, 5)


#-----------------------------------------------------------------------------#
# Extract top hits dosages from binarydosage 
#-----------------------------------------------------------------------------#
covariate_file <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_66485_GLM.rds")
dosages <- data.frame(readRDS("files/GetSNPValues_aspirin_age_ref_imp_sex_study_gxe_PC1-3_N_66485_chr5_out.rds")) %>% 
  rownames_to_column(var = 'vcfid')

dosages <- data.frame(readRDS("/home/rak/Dropbox/FIGI/Results/aspirin_190605/files/GetSNPValues_aspirin_age_ref_imp_sex_study_gxe_PC1-3_N_66485_chr5_out.rds")) %>% 
  rownames_to_column(var = 'vcfid')


# need to merge in dosage info to the EpiData
posthoc_df <- inner_join(covariate_file, dosages, by = 'vcfid')


# then apply gxe only to dosage columns, with the same covariates used for gxescan



# then compile information, create forest plot of betas CIs. 






















