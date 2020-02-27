#=============================================================================#
# Permutation testing for aspirin clumped finding... 
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())

exposure = 'aspirin'
covariates = c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3')
filename <- 'FIGI_v2.3_gxeset_aspirin_basic_covars_gxescan'

# set working directory
# important to do because figifs functions outputs into folders based on current position (e.g. figures)
setwd(paste0("~/Dropbox/FIGI/Results/", exposure))
gxe <- readRDS(paste0("~/data/results/", exposure, "/processed/", filename, "_results.rds"))
tmp <- do.call(rbind, lapply(list.files(paste0("~/data/results/", exposure, "/clump/"), full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))
gxe_chiSqGxE_ld <- gxe %>%
  filter(SNP %in% tmp$SNP)
rm(tmp)
gxe_twostep <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05)



# focus on bin #3 (has significant findings)
gxe_sig <- filter(gxe_twostep, grp == 3)

# need 'unclumped' set of SNPs for that bin.. 
table(gxe_sig$Chromosome)

unclump_wrap <- function(snp) {
  chr <- strsplit(snp, split = ":")[[1]][1] # get chromosome number
  tmp <- fread(paste0("~/data/results/aspirin/clump/FIGI_v2.3_gxeset_aspirin_basic_covars_gxescan_chr", chr, "_chiSqGxE_ldclump.clumped")) %>% 
    filter(SNP == snp)
  snp_list <- gsub("\\(1\\)", "", tmp$SP2)
  snp_list <- unlist(strsplit(snp_list, split = ","))
  snp_list
}

test <- sort(gxe_sig$SNP)
snps_in_ld <- unlist(lapply(test, unclump_wrap)) # make sure you also extract the tag SNPs and merge with this set
snps_in_ld <- snps_in_ld[snps_in_ld != "NONE"]
snps_in_ld_final <- sort(c(snps_in_ld, test))
snps_in_ld_final



#
# important - might have to redo above, using chiSqG clumping... because this has to do with bin assignemtn. 
#
#
 
#--------------------------------------#
# extract SNPs on HPC. then resume here
#
#--------------------------------------#

x1 <- readRDS("~/data/permutation_aspirin/FIGI_genotype_dosages_aspirin_permutation_1.rds")
x2 <- readRDS("~/data/permutation_aspirin/FIGI_genotype_dosages_aspirin_permutation_2.rds")
x3 <- readRDS("~/data/permutation_aspirin/FIGI_genotype_dosages_aspirin_permutation_3.rds")
x4 <- readRDS("~/data/permutation_aspirin/FIGI_genotype_dosages_aspirin_permutation_4.rds")
x5 <- readRDS("~/data/permutation_aspirin/FIGI_genotype_dosages_aspirin_permutation_5.rds")

x <- inner_join(x1, x2, 'vcfid')
x <- inner_join(x, x3, 'vcfid')
x <- inner_join(x, x4, 'vcfid')
x <- inner_join(x, x5, 'vcfid')

aspirin <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_aspirin_basic_covars_glm.rds")

perm_data <- inner_join(aspirin, x, 'vcfid')
colnames(perm_data) <- gsub("_.*$", "", colnames(perm_data))


snps_list_clump <- paste0("X", gsub(":", ".", test))
snps_list_noclump <- paste0("X", gsub(":", ".", snps_in_ld))

perm_data_clump <- perm_data[, ! names(perm_data) %in% snps_list_noclump]
perm_data_noclump <- perm_data[, ! names(perm_data) %in% snps_list_clump]

saveRDS(perm_data_clump, "~/data/permutation_aspirin/perm_data_clumped.rds")
saveRDS(perm_data, "~/data/permutation_aspirin/perm_data_not_clumped.rds")



model <- glm(outcome ~ 	X5.40238857.A.C*aspirin + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = perm_data, family = 'binomial')
summary(model)

# ok
model_base <- glm(outcome ~ X5.40238857.A.C+aspirin +nsaids+ age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = perm_data, family = 'binomial')
model_int  <- glm(outcome ~ X5.40238857.A.C*aspirin +nsaids+ age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = perm_data, family = 'binomial')
lrtest(model_base, model_int)


saveRDS(perm_data, "~/data/permutation_aspirin/perm_data.rds")


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



glm_func <- function(y) glm(outcome ~ y * aspirin + nsaids + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = perm_data, family = binomial)
glm_func_base <- function(y) glm(outcome ~ y + aspirin + nsaids + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = perm_data, family = binomial)
run_lrtest <- function(x,y) lrtest(x, y)


# run GLM + lrtest
basemodel <- map(perm_data[ , 10:19], function(x) glm_func_base(x))
intmodel <- map(perm_data[ , 10:19], function(x) glm_func(x))
results_gxe_glm <- mapply(run_lrtest, basemodel, intmodel, SIMPLIFY = F)

results_gxe_out <- do.call(rbind, results_gxe_glm) %>%
  tibble::rownames_to_column('SNP') %>%
  filter(!is.na(Chisq)) %>%
  mutate(SNP = gsub('.{2}$', '', SNP))


