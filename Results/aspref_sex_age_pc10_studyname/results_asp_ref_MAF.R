#========================================================#
# I want to calculate MAF on the GxE Top hits... 
# if you can also calculate by study..
#
#========================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(broom)
library(reshape)
rm(list = ls())

# results
results <- do.call(rbind, lapply(list.files(path = "./results/", full.names = T, pattern = "results_asp_ref_sex_age_pcs_studyname_N72438_"), fread, stringsAsFactors = F))
names(results)

results_GxE <- results %>%
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = F) %>% 
  mutate(CHR = as.numeric(CHR), 
         BP = as.numeric(BP), 
         P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, CHR, BP, Subjects, Cases, Reference, Alternate, betaGxE, P)

results_sig_GxE <- filter(results_GxE, P < 5e-8) # chromosomes 2, 3, 6, 7, 10, 16, 17, 18, 21, 22

results_GxE_filter <- results_GxE %>% 
  mutate(logP = -log10(P)) %>% 
  filter(!SNP %in% dupped$SNP) %>% arrange(BP)

# identify 'multiallelic' markers (or duplicates rows - maybe small gxescanR bug that captures rows twice)
dups <- results_GxE[duplicated(results_GxE$SNP), ]
dupped <- filter(results_GxE, SNP %in% dups$SNP)


# explore the top 10 GxE hits by studyname, etc. 
load("~/git/DATA/FIGI_EpiData_rdata/FIGI_Genotype_Epi_181120.RData")
pc30k <- fread("~/Dropbox/code/FIGI_Results/PrincipalComponents/FIGI_GxESet.eigenvec", skip = 1, 
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
  dplyr::select(vcfid, outcome, age_ref_imp, sex, studyname, paste0(rep("PC", 10), seq(1,10)), asp_ref, V2, platform) %>% 
  filter(complete.cases(.))

table(cov$studyname, cov$outcome)
exclude_studies <- c("ColoCare_1", "ESTHER_VERDI", "GALEON", "MOFFITT") # case-only studies.. 
cov <- filter(cov, !studyname %in% exclude_studies)

snps <- dplyr::select(cov, vcfid)
for(chr in c(2,3,6,7,10, 16, 17, 18, 21, 22)) {
  x <- as.data.frame(readRDS(paste0("/home/rak/Dropbox/code/FIGI_Results/aspref_sex_age_pc10_studyname/extract/GxEScanR_GxE_sig_loci_extract_chr", chr, ".rds"))) %>% rownames_to_column(var = 'vcfid')
  snps <- inner_join(snps, x, by = 'vcfid')
}
df <- inner_join(cov, snps, by = "vcfid")

# calculate MAFs
names(df)

# (alt allele freq, N = 72438)
df %>%summarise_at(c("2:4964187", "3:98981063",  "6:138017298", "7:99995536",  "10:3001065",  "16:79591072", "17:18528708", "18:3780526",  "21:34169318", "22:40651008"), function(x) mean(x))

af_studyname <- df %>% 
  group_by(studyname) %>% 
  summarise_at(c("2:4964187", "3:98981063",  "6:138017298", "7:99995536",  "10:3001065",  "16:79591072", "17:18528708", "18:3780526",  "21:34169318", "22:40651008"), function(x) mean(x))

af_studyname <- df %>% 
  group_by(platform) %>% 
  summarise_at(c("2:4964187", "3:98981063",  "6:138017298", "7:99995536",  "10:3001065",  "16:79591072", "17:18528708", "18:3780526",  "21:34169318", "22:40651008"), function(x) mean(x))

plot(af_studyname$`2:4964187`)
plot(af_studyname$`3:98981063`)
plot(af_studyname$`6:138017298`)
plot(af_studyname$`7:99995536`)
plot(af_studyname$`10:3001065`)
plot(af_studyname$`16:79591072`)
plot(af_studyname$`17:18528708`)
plot(af_studyname$`18:3780526`)
plot(af_studyname$`21:34169318`)
plot(af_studyname$`22:40651008`)

df %>%summarise_at(c("2:4964187", "3:98981063",  "6:138017298", "7:99995536",  "10:3001065",  "16:79591072", "17:18528708", "18:3780526",  "21:34169318", "22:40651008"), function(x) sum(x)/144876)








########################### WTF WTF WTF #####################

snps <- dplyr::select(cov, vcfid)
for(chr in c(2)) {
  x <- as.data.frame(readRDS(paste0("/home/rak/git/Dosage_Extract/GxEScanR_asp_ref_extract_chr", chr, ".rds"))) %>% rownames_to_column(var = 'vcfid')
  snps <- inner_join(snps, x, by = 'vcfid')
}
df <- inner_join(cov, snps, by = "vcfid")

# calculate MAFs
names(df)

# (alt allele freq, N = 72438)
df %>%summarise_at(c("2:4964187", "3:98981063",  "6:138017298", "7:99995536",  "10:3001065",  "16:79591072", "17:18528708", "18:3780526",  "21:34169318", "22:40651008"), function(x) mean(x))

bigwtf <- df %>% group_by(V2) %>%  summarise_at(c("2:96780986" , "2:109260147" ,"2:131698423", "2:172524900" ,"2:172558924", "2:172589295", "2:172600310" ,"2:172634504"), function(x) mean(x))

af_studyname <- df %>% 
  group_by(V2) %>% 
  summarise_at(c("2:4964187", "3:98981063",  "6:138017298", "7:99995536",  "10:3001065",  "16:79591072", "17:18528708", "18:3780526",  "21:34169318", "22:40651008"), function(x) mean(x))

plot(bigwtf$`2:96780986`)
plot(bigwtf$`2:109260147`)
plot(bigwtf$`2:131698423`)
plot(bigwtf$`2:172524900`)
plot(bigwtf$`2:172558924`)
plot(bigwtf$`2:172589295`)
plot(bigwtf$`2:172600310`)
plot(bigwtf$`2:172634504`)



x <- data.frame(table(df$V2, df$platform)) %>% 
  filter(Freq!=0)


###################################################
###################################################


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
