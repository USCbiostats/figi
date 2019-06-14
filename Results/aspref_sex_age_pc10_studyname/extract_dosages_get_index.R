#========================================================#
# before doing anything, we need to confirm that 
# 1) ids+dosages match between bdose and vcf files 
# 2) results match between gxescanR and glm
#
# Based on test scan of:
# outc ~ age+sex+pc+studyname+asp_ref
#
# Extract top 91 dosages from bdose and vcf files
#
# Ensure individuals have the same dosage values
# Ensure that sig results would be equally enormous using GLM
#========================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(broom)
library(lmtest)

rm(list = ls())
setwd("~/Dropbox/code/FIGI_Results/aspref_sex_age_pc10_studyname/")
source("~/Dropbox/code/FIGI_Results/aspref_sex_age_pc10_studyname/functions.R")

# original results
results <- do.call(rbind, lapply(list.files(full.names = T, pattern = "results_asp_ref_sex_age_pcs_studyname_N72438_"), fread, stringsAsFactors = F))
names(results)

#--------------------------------------------------------#
# Marginal G Results
# Get sig SNPs on genome-wide level (N = 91)

results_G <- results %>%
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = F) %>% 
  mutate(CHR = as.numeric(CHR), 
         BP = as.numeric(BP), 
         P = pchisq(chiSqG, df = 1, lower.tail = F)) %>% 
  dplyr::select(SNP, CHR, BP, Reference, Alternate, Subjects, Cases, betaG, chiSqG, P) %>% 
  filter(P < 5e-8)


# output text file to help extract markers from VCF files (vcftools)
# CHR \t POSITION
vcfout <- results_G %>% 
  dplyr::select(CHR, BP)
write.table(vcfout, file = "GxEScanR_asp_ref_top91_vcf_CHR_BP.txt", quote = F, row.names = F, col.names = F, sep = '\t')


# output R object to help extract markers from bdose files (BinaryDosage package)
# vector of CHR:BP
for(chr in 1:22) {
  snpsToGet <- results_G[which(results_G$CHR == chr), 'SNP']
  saveRDS(snpsToGet, file = paste0("GxEScanR_asp_ref_top91_bdose_SNP_chr", chr, ".rds"))
}

# x <- readRDS("~/Dropbox/working/GxEScanR_asp_ref_top91_bdose_SNP_chr15.rds")
# table(results_G$CHR)



##########################################################
##########################################################

# after extracting dosages from both BDOSE and VCFFILES...
# Check that dosages are identical
# if there are multiallelic markers, the check won't work, but let's make sure 
# dosages are the same for the ones that do ok. 

# read in covariate file used in GxEScanR run
cov <- readRDS("~/git/GxEScanR_CovFiles/FIGI_GxESet_asp_ref_sex_age_pcs_studyname_72438.rds")

# BDose Dosages
for(chr in c(6,8,15,18)) {
  x <- as.data.frame(readRDS(paste0("/home/rak/results/extract_tophits/binarydosage_asp_ref_age_sex_pc_studyname_N72438/GxEScanR_asp_ref_extract_chr", chr, ".rds"))) %>% 
    rownames_to_column(var = 'vcfid')
  cov <- inner_join(cov, x, by = 'vcfid')
}


# VCF Files
library(vcfR)
setwd("/home/rak/results/extract_tophits/vcftools_asp_ref_age_sex_pc_studyname_N72438")
x <- c(1:8,10:12,14,15,17:20)
x <- c(1,3:8,10,12,14,15,17:20) # remove 2, 11
x <- c(6, 8, 15, 18) # or do the 4 most sig chromosomes

# wrappper
vcf_wrap <- function(batch) {
  tmp <- lapply(x, function(y) vcfR2tidy(read.vcfR(paste0(batch,"_chr", y, "_asp_ref_top1262.recode.vcf"), verbose = F), single_frame = T, info_fields = c("CHROM", "POS"), format_fields = c("DS"))[[1]] %>% mutate(ID = paste(CHROM, POS, sep = ":")) )
  return(dcast(do.call(rbind, tmp), Indiv ~ ID, value.var = "gt_DS") %>% 
           dplyr::rename(vcfid = Indiv) %>% 
           dplyr::select(vcfid, everything()))
}


axiom_acs_aus_nf <- vcf_wrap('axiom_acs_aus_nf')
axiom_mecc_cfr_ky <- vcf_wrap('axiom_mecc_cfr_ky')
ccfr_1m_1mduo_reimpute <- vcf_wrap('ccfr_1m_1mduo_reimpute')
ccfr_omni <- vcf_wrap('ccfr_omni')
corect_oncoarray <- vcf_wrap('corect_oncoarray')
corect_oncoarray_nonEUR_reimpute <- vcf_wrap('corect_oncoarray_nonEUR_reimpute')
corsa_axiom <- vcf_wrap('corsa_axiom')
cytosnp_comb <- vcf_wrap('cytosnp_comb')
dachs3 <- vcf_wrap('dachs3')
initial_comb_datasets <- vcf_wrap('initial_comb_datasets')
mecc <- vcf_wrap('mecc')
newfoundland_omniquad <- vcf_wrap('newfoundland_omniquad')
omni_comb <- vcf_wrap('omni_comb')
omniexpress_exomechip <- vcf_wrap('omniexpress_exomechip')
oncoarray_to_usc <- vcf_wrap('oncoarray_to_usc')
plco_3 <- vcf_wrap('plco_3')
reach <- vcf_wrap('reach')

#ukbiobank <- vcf_wrap('ukbiobank') # not working.... 
# compare dosages without it, if good the likely biobank is good
# yeah.. biobank will probably be pretty mess. let's give it a try for a single chromosome

x <- read.vcfR("ukbiobank_chr8_asp_ref_top1262.recode.vcf") # look at 15th snp (interesting values)
vcfR2tidy(x) # error
y <- extract_gt_tidy(x) %>%  filter(Key ==83)

fff <- filter(vcfout, CHR == 8)

vcfy <- filter(cov, vcfid %in% y$Indiv) %>% 
  dplyr::select(vcfid, `8:128441535`)

ukbiobank_check <- inner_join(y, vcfy, by = c('Indiv' = 'vcfid')) %>% 
  mutate(gt_DS_round = round(gt_DS, 4)) 
all(ukbiobank_check$`8:128441535`== ukbiobank_check$gt_DS_round) # ukbiobank is rounded to the 4th digit, identical values

save.image(file = "vcfR_objects.RData")
load("vcfR_objects.RData")


# perform checks (without ukbiobank)
vcf <- rbind(axiom_acs_aus_nf, axiom_mecc_cfr_ky, ccfr_1m_1mduo_reimpute, ccfr_omni, corect_oncoarray, corect_oncoarray_nonEUR_reimpute, corsa_axiom, cytosnp_comb, dachs3, initial_comb_datasets, mecc, newfoundland_omniquad, omni_comb, omniexpress_exomechip, oncoarray_to_usc, plco_3, reach) %>% 
  filter(vcfid %in% cov$vcfid) %>% 
  arrange(vcfid)

bdose <- cov %>% 
  filter(vcfid %in% vcf$vcfid) %>% 
  dplyr::select(names(vcf)) %>% 
  arrange(vcfid)

all(vcf == bdose) # Perfect


# easy to display sameness...(in case you want to include in Rmd)
rm(list=setdiff(ls(), c("bdose", "vcf", "cov")))
save.image(file = "vcf_bdose_comparison_objects.RData")


##########################################################
##########################################################

# Run GLM on these markers...
# you can run it using studyname as factor, or using the same covariate file as gxescanr
# use the same covariate file (see below)


#------ GxE Set - asp_ref (N = 102,792 --> 72,438) ------
load("~/git/DATA/FIGI_EpiData_rdata/FIGI_Genotype_Epi_181120.RData")
pc30k <- fread("~/Dropbox/code/FIGI_PCA_IBD_Results/FIGI_GxESet.eigenvec", skip = 1, 
               col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))
cov <- Epi %>%
  filter(drop == 0 & gxe == 1) %>% 
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1),
         studyname = ifelse(studyname %in% c("NFCCR_1", "NFCCR_2"), "NFCCR", studyname),
         studyname = ifelse(studyname == "Colo2&3", "Colo23", studyname),
         asp_ref = ifelse(asp_ref == "Yes", 1, ifelse(asp_ref == "No", 0, NA))) %>% 
  dplyr::select(vcfid, outcome, age_ref_imp, sex, studyname, paste0(rep("PC", 10), seq(1,10)), asp_ref) %>% 
  filter(complete.cases(.))

table(cov$studyname, cov$outcome)
exclude_studies <- c("ColoCare_1", "ESTHER_VERDI", "GALEON", "MOFFITT") # case-only studies.. 

cov <- filter(cov, !studyname %in% exclude_studies)

# BDose Dosages
for(chr in 1:22) {
  x <- as.data.frame(readRDS(paste0("~/git/Dosage_Extract/GxEScanR_asp_ref_extract_chr", chr, ".rds"))) %>% 
    rownames_to_column(var = 'vcfid')
  cov <- inner_join(cov, x, by = 'vcfid')
}

# set studyname as a factor
cov <- cov %>% 
  mutate(studyname = as.factor(studyname))


#===============================================================#

# use the original one , to ensure there's no doubt
# read in covariate file used in GxEScanR run
cov <- readRDS("~/git/GxEScanR_CovFiles/FIGI_GxESet_asp_ref_sex_age_pcs_studyname_72438.rds")

# BDose Dosages - for this check, do a single chromosome
for(chr in c(8)) {
  x <- as.data.frame(readRDS(paste0("/home/rak/results/extract_tophits/binarydosage_asp_ref_age_sex_pc_studyname_N72438/GxEScanR_asp_ref_extract_chr", chr, ".rds"))) %>% 
    rownames_to_column(var = 'vcfid')
  cov <- inner_join(cov, x, by = 'vcfid')
}


glm_func_original <- function(y) glm(outcome ~ y + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CPSII_1 + MCCS_1 + NFCCR + MECC_2 + CCFR_3 + Kentucky + CCFR_1 + PPS3 + CRCGEN + MECC_3 + CCFR_4 + MEC_2 + ATBC + PPS4 + USC_HRT_CRC + MCCS_2 + DALS_2 + PLCO_2 + WHI_2 + DACHS_1 + Colo23 + MEC_1 + VITAL + DACHS_3 + DALS_1 + PLCO_1 + WHI_1 + MECC_1 + DACHS_2 + HPFS_3_AD + HPFS_1 + HPFS_2 + NHS_3_AD + NHS_2 + NHS_1 + PHS + NHS_4 + HPFS_4 + CLUEII + NHS_5_AD + PLCO_4_AD + EDRN + NCCCSI + NCCCSII + WHI_3 + CPSII_2 + SMS_AD + HPFS_5_AD + SELECT + PLCO_3 + REACH_AD, data = cov, family = binomial)


# results_glm_original <- map_dfr(cov[,67:99], function(x) glm_func_original(x) %>% tidy(), .id = "SNP")
# results_glm_original_filter <- results_glm_original %>%
#   filter(term == "y") %>%
#   separate(SNP, into = c("CHR", "BP"), remove = F) %>%
#   mutate(CHR = as.numeric(CHR),
#          BP = as.numeric(BP),
#          P = p.value)

# perform likelihood ratio test
# obtain base model
results_glm_original_base <- glm(outcome ~ asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CPSII_1 + MCCS_1 + NFCCR + MECC_2 + CCFR_3 + Kentucky + CCFR_1 + PPS3 + CRCGEN + MECC_3 + CCFR_4 + MEC_2 + ATBC + PPS4 + USC_HRT_CRC + MCCS_2 + DALS_2 + PLCO_2 + WHI_2 + DACHS_1 + Colo23 + MEC_1 + VITAL + DACHS_3 + DALS_1 + PLCO_1 + WHI_1 + MECC_1 + DACHS_2 + HPFS_3_AD + HPFS_1 + HPFS_2 + NHS_3_AD + NHS_2 + NHS_1 + PHS + NHS_4 + HPFS_4 + CLUEII + NHS_5_AD + PLCO_4_AD + EDRN + NCCCSI + NCCCSII + WHI_3 + CPSII_2 + SMS_AD + HPFS_5_AD + SELECT + PLCO_3 + REACH_AD, data = cov, family = binomial)

results_glm_original_base_tidy <- tidy(results_glm_original_base) # if you want to see estimates.. 

# run glm on each SNP, keep glm objects in list
results_glm_original_lr <- map(cov[,67:160], function(x) glm_func_original(x))

# run lrtest from lmtest package.
run_lrtest <- function(x) lrtest(results_glm_original_base, x) # library lmtest
results_glm_original_lrtest <- lapply(results_glm_original_lr, run_lrtest)

results_glm_original_lrtest_clean <- do.call(rbind, results_glm_original_lrtest) %>% 
  tibble::rownames_to_column('SNP') %>% 
  filter(!is.na(Chisq)) %>% 
  mutate(SNP = gsub('.{2}$', '', SNP))

z <- inner_join(results_G, results_glm_original_lrtest_clean, by = 'SNP') %>% 
  mutate(beta_gxescan = betaG, 
         P_gxescan = P, 
         P_glm = `Pr(>Chisq)`) %>% 
  dplyr::select(SNP, Subjects, Cases, beta_gxescan, P_gxescan, P_glm)

saveRDS(z, file = "results_G_GLM_chr8only.rds")

plot(log(z$P_gxescan), log(z$P_glm))
manhattan(z, p = 'P')




#========================================================#
# Replicate GxE results using GLM
#========================================================#
results_GxE <- results %>%
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = F) %>% 
  mutate(CHR = as.numeric(CHR), 
         BP = as.numeric(BP), 
         P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, CHR, BP, Subjects, Cases, Reference, Alternate, betaGxE, P)

results_GxE_filter <- results_GxE %>% 
  filter(P < 5e-08)


cov <- readRDS("~/git/GxEScanR_CovFiles/FIGI_GxESet_asp_ref_sex_age_pcs_studyname_72438.rds")
for(chr in c(2,3,6,7,10,16,17,18,21,22)) {
  x <- as.data.frame(readRDS(paste0("GxEScanR_GxE_sig_loci_extract_chr", chr, ".rds"))) %>% 
    rownames_to_column(var = 'vcfid')
  cov <- inner_join(cov, x, by = 'vcfid')
}


glm_func_original <- function(y) glm(outcome ~ y*asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CPSII_1 + MCCS_1 + NFCCR + MECC_2 + CCFR_3 + Kentucky + CCFR_1 + PPS3 + CRCGEN + MECC_3 + CCFR_4 + MEC_2 + ATBC + PPS4 + USC_HRT_CRC + MCCS_2 + DALS_2 + PLCO_2 + WHI_2 + DACHS_1 + Colo23 + MEC_1 + VITAL + DACHS_3 + DALS_1 + PLCO_1 + WHI_1 + MECC_1 + DACHS_2 + HPFS_3_AD + HPFS_1 + HPFS_2 + NHS_3_AD + NHS_2 + NHS_1 + PHS + NHS_4 + HPFS_4 + CLUEII + NHS_5_AD + PLCO_4_AD + EDRN + NCCCSI + NCCCSII + WHI_3 + CPSII_2 + SMS_AD + HPFS_5_AD + SELECT + PLCO_3 + REACH_AD, data = cov, family = binomial)

results_glm_original_base <- function(y) glm(outcome ~ y + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CPSII_1 + MCCS_1 + NFCCR + MECC_2 + CCFR_3 + Kentucky + CCFR_1 + PPS3 + CRCGEN + MECC_3 + CCFR_4 + MEC_2 + ATBC + PPS4 + USC_HRT_CRC + MCCS_2 + DALS_2 + PLCO_2 + WHI_2 + DACHS_1 + Colo23 + MEC_1 + VITAL + DACHS_3 + DALS_1 + PLCO_1 + WHI_1 + MECC_1 + DACHS_2 + HPFS_3_AD + HPFS_1 + HPFS_2 + NHS_3_AD + NHS_2 + NHS_1 + PHS + NHS_4 + HPFS_4 + CLUEII + NHS_5_AD + PLCO_4_AD + EDRN + NCCCSI + NCCCSII + WHI_3 + CPSII_2 + SMS_AD + HPFS_5_AD + SELECT + PLCO_3 + REACH_AD, data = cov, family = binomial)

# results_glm_original_base <- glm(outcome ~ y + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CPSII_1 + MCCS_1 + NFCCR + MECC_2 + CCFR_3 + Kentucky + CCFR_1 + PPS3 + CRCGEN + MECC_3 + CCFR_4 + MEC_2 + ATBC + PPS4 + USC_HRT_CRC + MCCS_2 + DALS_2 + PLCO_2 + WHI_2 + DACHS_1 + Colo23 + MEC_1 + VITAL + DACHS_3 + DALS_1 + PLCO_1 + WHI_1 + MECC_1 + DACHS_2 + HPFS_3_AD + HPFS_1 + HPFS_2 + NHS_3_AD + NHS_2 + NHS_1 + PHS + NHS_4 + HPFS_4 + CLUEII + NHS_5_AD + PLCO_4_AD + EDRN + NCCCSI + NCCCSII + WHI_3 + CPSII_2 + SMS_AD + HPFS_5_AD + SELECT + PLCO_3 + REACH_AD, data = cov, family = binomial)

# results_glm_original_base_tidy <- tidy(results_glm_original_base) # if you want to see estimates.. 

# run glm on each SNP, keep glm objects in list
results_glm_original_lrbase <- map(cov[,67:76], function(x) results_glm_original_base(x))
results_glm_original_lr <- map(cov[,67:76], function(x) glm_func_original(x))

# run lrtest from lmtest package
run_lrtest <- function(x,y) lrtest(x, y)
results_gxe <- mapply(run_lrtest, results_glm_original_lrbase, results_glm_original_lr, SIMPLIFY = F)

results_gxe_clean <- do.call(rbind, results_gxe) %>% 
  tibble::rownames_to_column('SNP') %>% 
  filter(!is.na(Chisq)) %>% 
  mutate(SNP = gsub('.{2}$', '', SNP))

z <- inner_join(results_GxE, results_gxe_clean, by = 'SNP') %>% 
mutate(betaGxE_gxescan = betaGxE, 
       P.GxE_gxescan = P, 
       P.GxE_glm = `Pr(>Chisq)`) %>% 
  dplyr::select(SNP, Subjects, Cases, betaGxE_gxescan, P.GxE_gxescan, P.GxE_glm)

saveRDS(z, file = "results_GxE_GLM.rds")



plot(log(z$P), log(z$`Pr(>Chisq)`))
manhattan(z, p = 'P')




