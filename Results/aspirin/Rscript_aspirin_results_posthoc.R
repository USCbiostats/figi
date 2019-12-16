#=============================================================================#
# NSAIDS - aspirin results
# 06/01/2019
# 09/20/2019
# 
# Follow-up analyses
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())

#-----------------------------------------------------------------------------#
# input significant result dosages
#-----------------------------------------------------------------------------#

# ... 
# ... 
# ...


# confirm results with GLM 
dat <- do.call(inner_join, lapply(list.files("~/data/Results/aspirin/dosage/", pattern = "GetSNPValues_aspirin_sex_age_pc3_studygxe_index_positions_chr", full.names = T), function(x) data.frame(readRDS(x)) %>% rownames_to_column("vcfid"))) %>% 
  merge(readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_66485_GLM.rds"))

glm_func <- function(y) glm(outcome ~ y * aspirin + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = dat, family = binomial)
glm_func_base <- function(y) glm(outcome ~ y + aspirin + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = dat, family = binomial)
run_lrtest <- function(x,y) lrtest(x, y)

# run GLM + lrtest
basemodel <- map(dat[ , 2:4], function(x) glm_func_base(x))
intmodel <- map(dat[,2:4], function(x) glm_func(x))
results_gxe_glm <- mapply(run_lrtest, basemodel, intmodel, SIMPLIFY = F)

results_gxe_out <- do.call(rbind, results_gxe_glm) %>%
  tibble::rownames_to_column('SNP') %>%
  filter(!is.na(Chisq)) %>%
  mutate(SNP = gsub('.{2}$', '', SNP))



# for X6.32560631, might be driven by UKB1.. 
dat <- do.call(inner_join, lapply(list.files("~/data/Results/aspirin/dosage/", pattern = "GetSNPValues_aspirin_sex_age_pc3_studygxe_index_positions_chr", full.names = T), function(x) data.frame(readRDS(x)) %>% rownames_to_column("vcfid"))) %>% 
  merge(readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_66485_GLM.rds"))

dat <- filter(dat, study_gxe != "UKB_1")

glm_func <- function(y) glm(outcome ~ y * aspirin + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = dat, family = binomial)
glm_func_base <- function(y) glm(outcome ~ y + aspirin + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = dat, family = binomial)
run_lrtest <- function(x,y) lrtest(x, y)

basemodel <- map(dat[ , 2:4], function(x) glm_func_base(x))
intmodel <- map(dat[,2:4], function(x) glm_func(x))
results_gxe_glm <- mapply(run_lrtest, basemodel, intmodel, SIMPLIFY = F)

results_gxe_out <- do.call(rbind, results_gxe_glm) %>%
  tibble::rownames_to_column('SNP') %>%
  filter(!is.na(Chisq)) %>%
  mutate(SNP = gsub('.{2}$', '', SNP))





#-----------------------------------------------------------------------------#
# Plot MAFs ----
# GENERALIZE THIS WITH FUNCTION
#-----------------------------------------------------------------------------#
posthoc_df_maf <- dat %>%
  group_by(study_gxe) %>%
  summarise_at(vars(X5.40252294:X6.12577203), function(x) 0.5 - abs( (sum(x) / (2*length(x))) - 0.5))


ggplot(posthoc_df_maf) +
  geom_point(aes(y = study_gxe, x = X5.40252294)) +
  theme_bw() +
  xlim(0,0.5)


ggplot(posthoc_df_maf) +
  geom_point(aes(y = study_gxe, x = X6.12577203)) +
  theme_bw() +
  xlim(0,0.5)



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
locuszoom_gxe <- gxe %>%
  mutate(`P-value` = pchisq(chiSqGxE, df = 1, lower.tail = F),
         MarkerName = paste0("chr", SNP)) %>%
  dplyr::select(MarkerName, `P-value`)
write.table(locuszoom_gxe, file = "/media/work/tmp/LocusZoom_GxE_aspirin_sex_age_pc3_studygxe_66485.txt", quote = F, row.names = F, sep = "\t")


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






















