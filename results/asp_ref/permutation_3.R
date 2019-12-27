#=============================================================================#
# let's see if we can adapt this to gxe
# it's slightly different - you're not creating a distribution of a value like correlation coef
# instead you're counting number of significant hits @ some pre-determined cutoff value
#
# let's start with aspirin, since the significant finding is in bin #2
#=============================================================================#
library(tidyverse)
library(data.table)
library(figifs)
library(lmtest)


aspirin <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_72269_GLM.rds")


dosages <- readRDS("~/permutation/input/aspirin_bin2_hit_ld_snplist_index_TMP1.rds") %>% 
  filter(vcfid %in% aspirin$vcfid)

perm_data <- inner_join(aspirin, dosages, by = 'vcfid') %>% 
  dplyr::select(outcome, age_ref_imp, sex, study_gxe, PC1, PC2, PC3, aspirin, colnames(dosages))

saveRDS(perm_data, file = "~/permutation/input/perm_data2.rds", version = 2)

#-----------------------------------------------------------------------------#
# actual permutation?
#-----------------------------------------------------------------------------#

# each permutation consists of fitting models for all 10 SNPs, then counting the number of SNPs that are 'significant'. 

# pseudo code

# for each permutation - 
#   apply columns glm with same covariates
#   count if snps p value < sig number 1
#   append number to vector 1
#   count if snps p value < sig number 2
#   append number to vector 2
#   
# do this for each snp in each iteration (testing if its smaller than sig levels)
# so.. internal loop counter




# lrtest for a single chromosome



glm_func <- function(y) glm(outcome ~ y * aspirin + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = perm_data, family = binomial)
glm_func_base <- function(y) glm(outcome ~ y + aspirin + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = perm_data, family = binomial)
run_lrtest <- function(x,y) lrtest(x, y)

# run GLM + lrtest
basemodel <- map(perm_data[ , 10:19], function(x) glm_func_base(x))
intmodel <- map(perm_data[ , 10:19], function(x) glm_func(x))
results_gxe_glm <- mapply(run_lrtest, basemodel, intmodel, SIMPLIFY = F)

results_gxe_out <- do.call(rbind, results_gxe_glm) %>%
  tibble::rownames_to_column('SNP') %>%
  filter(!is.na(Chisq)) %>%
  mutate(SNP = gsub('.{2}$', '', SNP))



#-----------------------------------------------------------------------------#
# gather up results, create counts!
#-----------------------------------------------------------------------------#
filelist <- list.files("~/permutation/", pattern = "permutation_", full.names = T)

get_counts_binp <- function(x) {
  sum(readRDS(x)[, c("Pr(>Chisq)")] < 0.00125)
}


filelist <- list.files("~/permutation/", pattern = "permutation_", full.names = T)
out <- do.call(c, lapply(filelist, get_counts_binp))
sum(out)
sum(out)/length(filelist)
hist(out)


get_counts_binp <- function(x) {
  sum(readRDS(x)[, c("Pr(>Chisq)")] < 0.0001197766)
}
out <- do.call(c, lapply(filelist, get_counts_binp))

hist(out)

