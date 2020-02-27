
# first, make sure that you're outputting E|G statistics for locuszoom
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())
gxe <- readRDS("~/data/results/asp_ref/processed/FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan_results.rds")

# output E|G pvalues
locuszoom <- gxe %>%
  dplyr::mutate(`P-value` = pchisq(chiSqControl, df = 1, lower.tail = F),
                MarkerName = paste0("chr", Chromosome, ":", Location)) %>% 
  dplyr::select(MarkerName, `P-value`)

write.table(locuszoom, file = paste0("/media/work/tmp/", "FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_gxescan", "_chiSqControl_locuszoom.txt"), quote = F, row.names = F, sep = "\t")

# output E|G pvalues
locuszoom <- gxe %>%
  dplyr::mutate(`P-value` = pchisq(chiSqGE, df = 1, lower.tail = F),
                MarkerName = paste0("chr", Chromosome, ":", Location)) %>% 
  dplyr::select(MarkerName, `P-value`)

write.table(locuszoom, file = paste0("/media/work/tmp/", "FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_gxescan", "_chiSqGE_locuszoom.txt"), quote = F, row.names = F, sep = "\t")

# output 3df pvalues
locuszoom <- gxe %>%
  dplyr::mutate(`P-value` = pchisq(chiSq3df, df = 3, lower.tail = F),
                MarkerName = paste0("chr", Chromosome, ":", Location)) %>% 
  dplyr::select(MarkerName, `P-value`)

write.table(locuszoom, file = paste0("/media/work/tmp/", "FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_gxescan", "_chiSq3df_locuszoom.txt"), quote = F, row.names = F, sep = "\t")




# # original code to output GxE statistics
# locuszoom <- gxe %>%
#   dplyr::mutate(`P-value` = pchisq(chiSqGxE, df = 1, lower.tail = F),
#                 MarkerName = paste0("chr", Chromosome, ":", Location)) %>% 
#   dplyr::select(MarkerName, `P-value`)
# 
# write.table(locuszoom, file = paste0("/media/work/tmp/", filename, "_chiSqGxE_locuszoom.txt"), quote = F, row.names = F, sep = "\t")
# 
# # original code to output G statistics 
# locuszoom <- gxe %>%
#   dplyr::mutate(`P-value` = pchisq(chiSqG, df = 1, lower.tail = F),
#                 MarkerName = paste0("chr", Chromosome, ":", Location)) %>% 
#   dplyr::select(MarkerName, `P-value`)
# 
# write.table(locuszoom, file = paste0("/media/work/tmp/", filename, "_chiSqG_locuszoom.txt"), quote = F, row.names = F, sep = "\t")




# imputation quality and minor allele frequency
x <- readRDS("~/data/HRC_InfoFile_Merged/HRC_InfoFile_Chr6.rds") %>% 
  filter(ID == "6:12577203:T:C")

xplot <- mutate_at(x, vars(matches("Rsq")), as.numeric) %>% 
  pivot_longer( cols = contains("Rsq"), names_to = "study_gxe", values_to = "Rsq") %>% 
  mutate(study_gxe = gsub("_Rsq", "", study_gxe))

# names(xplot) <- c("ID","axiom_acs_aus_nf_ALT_Frq",       "axiom_mecc_cfr_ky_ALT_Frq"      ,"ccfr_1m_1mduo_reimpute_ALT_Frq", "ccfr_omni_ALT_Frq"  ,            "corect_oncoarray_ALT_Frq"   ,    "corsa_axiom_ALT_Frq"      ,      "cytosnp_comb_ALT_Frq" ,"initial_comb_datasets_ALT_Frq"  ,"mecc_ALT_Frq"   ,               "newfoundland_omniquad_ALT_Frq" , "omni_comb_ALT_Frq"     , "omniexpress_exomechip_ALT_Frq" ,"oncoarray_to_usc_ALT_Frq"    ,   "plco_3_ALT_Frq"         ,        "reach_ALT_Frq"        ,          "ukbiobank_ALT_Frq"        ,      "study_gxe"       ,               "Rsq")
# 


p<-ggplot(xplot, aes(Rsq, study_gxe)) + 
  geom_point(size = 3, color = 'red') + 
  theme_bw() + 
  theme(text=element_text(size = 24)) + 
  xlim(0,1)
p
ggsave("~/Dropbox/plot_6_12577203_rsq.png", plot = p)




xplot <- mutate_at(x, vars(matches("Frq")), as.numeric) %>% 
  pivot_longer( cols = contains("Frq"), names_to = "study_gxe", values_to = "MAF") %>% 
  # mutate(MAF = 1-MAF) %>% 
  mutate(study_gxe = gsub("_ALT_Frq", "", study_gxe))

p<-ggplot(xplot, aes(MAF, study_gxe)) + 
  geom_point(size = 3, color = 'blue') + 
  theme_bw() + 
  theme(text=element_text(size = 24)) + 
  xlim(0,0.5)
p
ggsave("~/Dropbox/plot_6_12577203_maf.png", plot = p)







# imputation quality and minor allele frequency
x <- readRDS("~/data/HRC_InfoFile_Merged/HRC_InfoFile_Chr6.rds") %>% 
  filter(ID == "6:32560631:C:T")

xplot <- mutate_at(x, vars(matches("Rsq")), as.numeric) %>% 
  pivot_longer( cols = contains("Rsq"), names_to = "study_gxe", values_to = "Rsq") %>% 
  mutate(study_gxe = gsub("_Rsq", "", study_gxe))

p<-ggplot(xplot, aes(Rsq, study_gxe)) + 
  geom_point(size = 3, color = 'red') + 
  theme_bw() + 
  theme(text=element_text(size = 24)) + 
  xlim(0,1)
p
ggsave("~/Dropbox/plot_6_32560631_rsq.png", plot = p)




xplot <- mutate_at(x, vars(matches("Frq")), as.numeric) %>% 
  pivot_longer( cols = contains("Frq"), names_to = "study_gxe", values_to = "MAF") %>% 
  # mutate(MAF = 1-MAF) %>%
  mutate(study_gxe = gsub("_ALT_Frq", "", study_gxe))

p<-ggplot(xplot, aes(MAF, study_gxe)) + 
  geom_point(size = 3, color = 'blue') + 
  theme_bw() + 
  theme(text=element_text(size = 24)) + 
  xlim(0,0.5)
p
ggsave("~/Dropbox/plot_6_32560631_maf.png", plot = p)







# let's explore 6_32560631 in terms of UKB and the rest of the cohorts

dat <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_asp_ref_basic_covars_glm.rds")
geno <- readRDS("~/FIGI_genotype_dosages_asp_ref_gxe.rds")

dat <- inner_join(dat, geno, 'vcfid')


dat_ukb <- filter(dat, study_gxe == "UKB_1")
dat_rest <- filter(dat, study_gxe != "UKB_1")

glm_all <- glm(outcome ~ X6.32560631.C.T * asp_ref + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = dat, family = 'binomial')
summary(glm_all)


glm_rest <- glm(outcome ~ X6.32560631.C.T * asp_ref + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = dat_rest, family = 'binomial')
summary(glm_rest)


glm_ukb <- glm(outcome ~ X6.32560631.C.T * asp_ref + age_ref_imp + sex              + PC1 + PC2 + PC3, data = dat_ukb, family = 'binomial')
summary(glm_ukb)

# main effects of the gene

glm_rest <- glm(outcome ~ X6.32560631.C.T +  asp_ref + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = dat_rest, family = 'binomial')
summary(glm_rest)


glm_ukb <- glm(outcome ~ X6.32560631.C.T +  asp_ref + age_ref_imp + sex              + PC1 + PC2 + PC3, data = dat_ukb, family = 'binomial')
summary(glm_ukb)



# main effects of the exposure

glm_rest <- glm(outcome ~  asp_ref + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = dat_rest, family = 'binomial')
summary(glm_rest)


glm_ukb <- glm(outcome ~  asp_ref + age_ref_imp + sex              + PC1 + PC2 + PC3, data = dat_ukb, family = 'binomial')
summary(glm_ukb)



# explore G|E associations by study
glm_all <- glm(asp_ref ~ X6.32560631.C.T + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = dat)
summary(glm_all)

glm_rest <- glm(asp_ref ~ X6.32560631.C.T + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = dat_rest)
summary(glm_rest)

glm_ukb <- glm(asp_ref ~ X6.32560631.C.T + age_ref_imp + sex + PC1 + PC2 + PC3, data = dat_ukb)
summary(glm_ukb)


cor(dat_rest$X6.32560631.C.T, dat_rest$asp_ref)
cor(dat_ukb$X6.32560631.C.T, dat_ukb$asp_ref)



# distribution of asp_ref by ukb, not-ukb
tmp <- dat %>% 
  mutate(ukb = ifelse(study_gxe == "UKB_1", "UKB", "non UKB"))

table(tmp$asp_ref, tmp$ukb)

table(round(tmp$X6.32560631.C.T, 1))



# later on - meta-analysis
posthoc_df <- do.call(inner_join, lapply(list.files("~/data/Results/aspirin/dosage/", pattern = "GetSNPValues_aspirin_sex_age_pc3_studygxe_index_positions_chr", full.names = T), function(x) data.frame(readRDS(x)))) %>% 
  merge(readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_72269_GLM.rds"))

study_info <- readRDS("~/data/Annotations/FIGI_studygxe_info.rds")
gxe_counts <- get_counts_outcome_by_group(posthoc_df, outcome, X5.40252294, study_gxe)
gxe_glm <- get_estimates_gxe_by_group(posthoc_df, outcome, aspirin, study_gxe, X5.40252294, age_ref_imp, sex)
gxe_meta <- dplyr::bind_cols(gxe_counts, gxe_glm) %>%
  inner_join(., study_info, by = "study_gxe")

results_meta <- meta::metagen(estimate,
                              std.error,
                              data=gxe_meta,
                              studlab=paste(study_gxe),
                              comb.fixed = FALSE,
                              comb.random = TRUE,
                              method.tau = "SJ",
                              hakn = TRUE,
                              prediction=TRUE,
                              sm="OR",
                              byvar = study_design)

meta::forest(results_meta,
             layout = "JAMA",
             text.predict = "95% CI",
             col.predict = "black",
             leftcols = c("studlab", "Control", "Case", "N", "effect", "ci", "w.random"),
             digits.addcols=0)

meta::funnel(results_meta, sm="OR", studlab = T)

