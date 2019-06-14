#========================================================#
# explore effects of different PC adjustments
#
# for convenience, use chr8 only
#========================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(broom)
rm(list = ls())

# results
results <- do.call(rbind, lapply(list.files(full.names = T, pattern = "results_asp_ref_sex_age_pcs_studyname_N72438_"), fread, stringsAsFactors = F))
names(results)

results_G <- results %>%
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = F) %>% 
  mutate(CHR = as.numeric(CHR), 
         BP = as.numeric(BP), 
         P = pchisq(chiSqG, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, CHR, BP, Subjects, Cases, Reference, Alternate, betaG, P)

# chromosome 8..
results_G_filter <- results_G %>% 
  mutate(logP = -log10(P)) %>% 
  filter(CHR == 8,
         logP >= 27.49965) # get top 10 for this comparison

results_G_filter <- results_G %>% 
  mutate(logP = -log10(P)) %>% 
  filter(CHR == 8) %>%  arrange(desc(logP)) # get top 10 for this comparison

write.table(results_G_filter, file = "~/Dropbox/test.txt", quote = F, row.names = F)

#--------------------------------------------------------#
# replicate analysis using GLM on select top hits on chromosome 8 (arbitrary choice, chromosome with very high significance hit)
cov <- readRDS("~/git/GxEScanR_CovFiles/FIGI_GxESet_asp_ref_sex_age_pcs20_studyname_72438.rds")
snps <- dplyr::select(cov, vcfid)
for(chr in c(8)) {
  x <- as.data.frame(readRDS(paste0("/home/rak/results/extract_tophits/binarydosage_asp_ref_age_sex_pc_studyname_N72438/GxEScanR_asp_ref_extract_chr", chr, ".rds"))) %>% rownames_to_column(var = 'vcfid')
  snps <- inner_join(snps, x, by = 'vcfid')
}
snps <- dplyr::select(snps, one_of(c('vcfid', results_G_filter$SNP)))
df <- inner_join(cov, snps, by = "vcfid")

# for comparison of different adjustment models (PC exploration)
glm_function_0 <- function(y) glm(outcome ~ y + asp_ref + age_ref_imp + sex + CPSII_1 + MCCS_1 + NFCCR + MECC_2 + CCFR_3 + Kentucky + CCFR_1 + PPS3 + CRCGEN + MECC_3 + CCFR_4 + MEC_2 + ATBC + PPS4 + USC_HRT_CRC + MCCS_2 + DALS_2 + PLCO_2 + WHI_2 + DACHS_1 + Colo23 + MEC_1 + VITAL + DACHS_3 + DALS_1 + PLCO_1 + WHI_1 + MECC_1 + DACHS_2 + HPFS_3_AD + HPFS_1 + HPFS_2 + NHS_3_AD + NHS_2 + NHS_1 + PHS + NHS_4 + HPFS_4 + CLUEII + NHS_5_AD + PLCO_4_AD + EDRN + NCCCSI + NCCCSII + WHI_3 + CPSII_2 + SMS_AD + HPFS_5_AD + SELECT + PLCO_3 + REACH_AD, data = df, family = binomial)
glm_function_4 <- function(y) glm(outcome ~ y + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + CPSII_1 + MCCS_1 + NFCCR + MECC_2 + CCFR_3 + Kentucky + CCFR_1 + PPS3 + CRCGEN + MECC_3 + CCFR_4 + MEC_2 + ATBC + PPS4 + USC_HRT_CRC + MCCS_2 + DALS_2 + PLCO_2 + WHI_2 + DACHS_1 + Colo23 + MEC_1 + VITAL + DACHS_3 + DALS_1 + PLCO_1 + WHI_1 + MECC_1 + DACHS_2 + HPFS_3_AD + HPFS_1 + HPFS_2 + NHS_3_AD + NHS_2 + NHS_1 + PHS + NHS_4 + HPFS_4 + CLUEII + NHS_5_AD + PLCO_4_AD + EDRN + NCCCSI + NCCCSII + WHI_3 + CPSII_2 + SMS_AD + HPFS_5_AD + SELECT + PLCO_3 + REACH_AD, data = df, family = binomial)
glm_function_10 <- function(y) glm(outcome ~ y + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CPSII_1 + MCCS_1 + NFCCR + MECC_2 + CCFR_3 + Kentucky + CCFR_1 + PPS3 + CRCGEN + MECC_3 + CCFR_4 + MEC_2 + ATBC + PPS4 + USC_HRT_CRC + MCCS_2 + DALS_2 + PLCO_2 + WHI_2 + DACHS_1 + Colo23 + MEC_1 + VITAL + DACHS_3 + DALS_1 + PLCO_1 + WHI_1 + MECC_1 + DACHS_2 + HPFS_3_AD + HPFS_1 + HPFS_2 + NHS_3_AD + NHS_2 + NHS_1 + PHS + NHS_4 + HPFS_4 + CLUEII + NHS_5_AD + PLCO_4_AD + EDRN + NCCCSI + NCCCSII + WHI_3 + CPSII_2 + SMS_AD + HPFS_5_AD + SELECT + PLCO_3 + REACH_AD, data = df, family = binomial)
glm_function_20 <- function(y) glm(outcome ~ y + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + CPSII_1 + MCCS_1 + NFCCR + MECC_2 + CCFR_3 + Kentucky + CCFR_1 + PPS3 + CRCGEN + MECC_3 + CCFR_4 + MEC_2 + ATBC + PPS4 + USC_HRT_CRC + MCCS_2 + DALS_2 + PLCO_2 + WHI_2 + DACHS_1 + Colo23 + MEC_1 + VITAL + DACHS_3 + DALS_1 + PLCO_1 + WHI_1 + MECC_1 + DACHS_2 + HPFS_3_AD + HPFS_1 + HPFS_2 + NHS_3_AD + NHS_2 + NHS_1 + PHS + NHS_4 + HPFS_4 + CLUEII + NHS_5_AD + PLCO_4_AD + EDRN + NCCCSI + NCCCSII + WHI_3 + CPSII_2 + SMS_AD + HPFS_5_AD + SELECT + PLCO_3 + REACH_AD, data = df, family = binomial)

# run GLM on each SNP
results_dosage_pc0 <- map(df[,77:86], function(x) glm_function_0(x))
results_dosage_pc4 <- map(df[,77:86], function(x) glm_function_4(x))
results_dosage_pc10 <- map(df[,77:86], function(x) glm_function_10(x))
results_dosage_pc20 <- map(df[,77:86], function(x) glm_function_20(x))



results_dosage_pc0_tidy <- do.call(rbind, lapply(results_dosage_pc0, tidy)) %>% 
  rownames_to_column('SNP') %>%
  mutate(SNP = gsub('.{2}$', '', SNP),
         est_pc0 = estimate) %>% 
  filter(term == "y")
results_dosage_pc4_tidy <- do.call(rbind, lapply(results_dosage_pc4, tidy)) %>% 
  rownames_to_column('SNP') %>%
  mutate(SNP = gsub('.{2}$', '', SNP),
         est_pc4 = estimate) %>% 
  filter(term == "y")
results_dosage_pc10_tidy <- do.call(rbind, lapply(results_dosage_pc10, tidy)) %>% 
  rownames_to_column('SNP') %>%
  mutate(SNP = gsub('.{2}$', '', SNP),
         est_pc10 = estimate) %>% 
  filter(term == "y")
results_dosage_pc20_tidy <- do.call(rbind, lapply(results_dosage_pc20, tidy)) %>% 
  rownames_to_column('SNP') %>%
  mutate(SNP = gsub('.{2}$', '', SNP),
         est_pc20 = estimate) %>% 
  filter(term == "y")






# base models for likelihood ratio
results_glm_original_base_pc0 <- glm(outcome ~ asp_ref + age_ref_imp + sex + CPSII_1 + MCCS_1 + NFCCR + MECC_2 + CCFR_3 + Kentucky + CCFR_1 + PPS3 + CRCGEN + MECC_3 + CCFR_4 + MEC_2 + ATBC + PPS4 + USC_HRT_CRC + MCCS_2 + DALS_2 + PLCO_2 + WHI_2 + DACHS_1 + Colo23 + MEC_1 + VITAL + DACHS_3 + DALS_1 + PLCO_1 + WHI_1 + MECC_1 + DACHS_2 + HPFS_3_AD + HPFS_1 + HPFS_2 + NHS_3_AD + NHS_2 + NHS_1 + PHS + NHS_4 + HPFS_4 + CLUEII + NHS_5_AD + PLCO_4_AD + EDRN + NCCCSI + NCCCSII + WHI_3 + CPSII_2 + SMS_AD + HPFS_5_AD + SELECT + PLCO_3 + REACH_AD, data = df, family = binomial)
results_glm_original_base_pc4 <- glm(outcome ~ asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + CPSII_1 + MCCS_1 + NFCCR + MECC_2 + CCFR_3 + Kentucky + CCFR_1 + PPS3 + CRCGEN + MECC_3 + CCFR_4 + MEC_2 + ATBC + PPS4 + USC_HRT_CRC + MCCS_2 + DALS_2 + PLCO_2 + WHI_2 + DACHS_1 + Colo23 + MEC_1 + VITAL + DACHS_3 + DALS_1 + PLCO_1 + WHI_1 + MECC_1 + DACHS_2 + HPFS_3_AD + HPFS_1 + HPFS_2 + NHS_3_AD + NHS_2 + NHS_1 + PHS + NHS_4 + HPFS_4 + CLUEII + NHS_5_AD + PLCO_4_AD + EDRN + NCCCSI + NCCCSII + WHI_3 + CPSII_2 + SMS_AD + HPFS_5_AD + SELECT + PLCO_3 + REACH_AD, data = df, family = binomial)
results_glm_original_base_pc10 <- glm(outcome ~ asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CPSII_1 + MCCS_1 + NFCCR + MECC_2 + CCFR_3 + Kentucky + CCFR_1 + PPS3 + CRCGEN + MECC_3 + CCFR_4 + MEC_2 + ATBC + PPS4 + USC_HRT_CRC + MCCS_2 + DALS_2 + PLCO_2 + WHI_2 + DACHS_1 + Colo23 + MEC_1 + VITAL + DACHS_3 + DALS_1 + PLCO_1 + WHI_1 + MECC_1 + DACHS_2 + HPFS_3_AD + HPFS_1 + HPFS_2 + NHS_3_AD + NHS_2 + NHS_1 + PHS + NHS_4 + HPFS_4 + CLUEII + NHS_5_AD + PLCO_4_AD + EDRN + NCCCSI + NCCCSII + WHI_3 + CPSII_2 + SMS_AD + HPFS_5_AD + SELECT + PLCO_3 + REACH_AD, data = df, family = binomial)
results_glm_original_base_pc20 <- glm(outcome ~ asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + CPSII_1 + MCCS_1 + NFCCR + MECC_2 + CCFR_3 + Kentucky + CCFR_1 + PPS3 + CRCGEN + MECC_3 + CCFR_4 + MEC_2 + ATBC + PPS4 + USC_HRT_CRC + MCCS_2 + DALS_2 + PLCO_2 + WHI_2 + DACHS_1 + Colo23 + MEC_1 + VITAL + DACHS_3 + DALS_1 + PLCO_1 + WHI_1 + MECC_1 + DACHS_2 + HPFS_3_AD + HPFS_1 + HPFS_2 + NHS_3_AD + NHS_2 + NHS_1 + PHS + NHS_4 + HPFS_4 + CLUEII + NHS_5_AD + PLCO_4_AD + EDRN + NCCCSI + NCCCSII + WHI_3 + CPSII_2 + SMS_AD + HPFS_5_AD + SELECT + PLCO_3 + REACH_AD, data = df, family = binomial)






# get lrtest P VALUES

# run lrtest from lmtest package
# no PC adjustment
run_lrtest <- function(x) lrtest(results_glm_original_base_pc0, x) # library lmtest
results_lrtest <- lapply(results_dosage_pc0, run_lrtest)
results_lrtest_clean_pc0 <- do.call(rbind, results_lrtest) %>% 
  tibble::rownames_to_column('SNP') %>% 
  filter(!is.na(Chisq)) %>% 
  mutate(SNP = gsub('.{2}$', '', SNP),
         P_pc0 = `Pr(>Chisq)`,
         logP_pc0 = -log10(`Pr(>Chisq)`))


run_lrtest <- function(x) lrtest(results_glm_original_base_pc4, x) # library lmtest
results_lrtest <- lapply(results_dosage_pc4, run_lrtest)
results_lrtest_clean_pc4 <- do.call(rbind, results_lrtest) %>% 
  tibble::rownames_to_column('SNP') %>% 
  filter(!is.na(Chisq)) %>% 
  mutate(SNP = gsub('.{2}$', '', SNP),
         P_pc4 = `Pr(>Chisq)`,
         logP_pc4 = -log10(`Pr(>Chisq)`))

run_lrtest <- function(x) lrtest(results_glm_original_base_pc10, x) # library lmtest
results_lrtest <- lapply(results_dosage_pc10, run_lrtest)
results_lrtest_clean_pc10 <- do.call(rbind, results_lrtest) %>% 
  tibble::rownames_to_column('SNP') %>% 
  filter(!is.na(Chisq)) %>% 
  mutate(SNP = gsub('.{2}$', '', SNP),
         P_pc10 = `Pr(>Chisq)`,
         logP_pc10 = -log10(`Pr(>Chisq)`))


run_lrtest <- function(x) lrtest(results_glm_original_base_pc20, x) # library lmtest
results_lrtest <- lapply(results_dosage_pc20, run_lrtest)
results_lrtest_clean_pc20 <- do.call(rbind, results_lrtest) %>% 
  tibble::rownames_to_column('SNP') %>% 
  filter(!is.na(Chisq)) %>% 
  mutate(SNP = gsub('.{2}$', '', SNP),
         P_pc20 = `Pr(>Chisq)`,
         logP_pc20 = -log10(`Pr(>Chisq)`))


z <- cbind(results_G_filter[, c('SNP', 'CHR', 'BP', 'Subjects', 'Cases')], 
               results_dosage_pc0_tidy[, c('est_pc0')], 
               results_lrtest_clean_pc0[, c('P_pc0', 'logP_pc0')], 
               results_dosage_pc4_tidy[, c('est_pc4')], 
               results_lrtest_clean_pc4[, c('P_pc4', 'logP_pc4')],
               results_dosage_pc10_tidy[, c('est_pc10')], 
               results_lrtest_clean_pc10[, c('P_pc10', 'logP_pc10')],
               results_dosage_pc20_tidy[, c('est_pc20')], 
               results_lrtest_clean_pc20[, c('P_pc20', 'logP_pc20')]) 

z_ests <- cbind(results_G_filter[, c('SNP', 'CHR', 'BP', 'Subjects', 'Cases')], 
           results_dosage_pc0_tidy[, c('est_pc0')], 
           results_dosage_pc4_tidy[, c('est_pc4')], 
           results_dosage_pc10_tidy[, c('est_pc10')], 
           results_dosage_pc20_tidy[, c('est_pc20')])
saveRDS(z_ests, file = "~/z_ests.rds")


z_pvals <- cbind(results_G_filter[, c('SNP', 'CHR', 'BP', 'Subjects', 'Cases')],
           results_lrtest_clean_pc0[, c('P_pc0', 'logP_pc0')], 
           results_lrtest_clean_pc4[, c('P_pc4', 'logP_pc4')],
           results_lrtest_clean_pc10[, c('P_pc10', 'logP_pc10')],
           results_lrtest_clean_pc20[, c('P_pc20', 'logP_pc20')]) 

z_pvals <- cbind(results_G_filter[, c('SNP', 'CHR', 'BP', 'Subjects', 'Cases')],
                 results_lrtest_clean_pc0[, c('logP_pc0')], 
                 results_lrtest_clean_pc4[, c('logP_pc4')],
                 results_lrtest_clean_pc10[, c('logP_pc10')],
                 results_lrtest_clean_pc20[, c('logP_pc20')]) 
saveRDS(z_pvals, file = "~/z_pvals.rds")


par(mfrow = c(2, 2))
plot(z$logP_pc10, z$logP_pc0, main = "logP noPC vs PC1-10")
plot(z$logP_pc10, z$logP_pc4, main = "logP PC1-4 vs PC1-10")
plot(z$logP_pc10, z$logP_pc20, main = "logP PC1-20 vs PC1-10")

par(mfrow = c(2, 2))
plot(z$BP, z$logP_pc0, main = "G main effects, Chr 8, no PC adj")
plot(z$BP, z$logP_pc4, main = "G main effects, Chr 8, PC 1-4")
plot(z$BP, z$logP_pc10, main = "G main effects, Chr 8, PC 1-10")
plot(z$BP, z$logP_pc20, main = "G main effects, Chr 8, PC 1-20")


  # mutate(beta_gxescan = betaG, 
  #        P_gxescan = P, 
  #        P_glm = `Pr(>Chisq)`) %>% 
  # dplyr::select(SNP, Subjects, Cases, beta_gxescan, P_gxescan, P_glm)

saveRDS(z, file = "results_G_GLM_chr8only.rds")

plot(log(z$P_gxescan), log(z$P_glm))
manhattan(z, p = 'P')