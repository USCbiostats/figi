# for figi results reports, need a bunch of additional things
# let's work on those here


# LD structure
# I think it's informative to use locuszoom plots for this, based on the E|G results among controls (or all samples, not sure which to use)

# first, make sure that you're outputting E|G statistics for locuszoom
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())
gxe <- readRDS("~/data/results/alcoholc_moderate/processed/FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_gxescan_results.rds")

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
x <- readRDS("~/data/HRC_InfoFile_Merged/HRC_InfoFile_Chr17.rds") %>% 
  filter(ID == "17:27861663:G:A")

xplot <- mutate_at(x, vars(matches("Rsq")), as.numeric) %>% 
  pivot_longer( cols = contains("Rsq"), names_to = "study_gxe", values_to = "Rsq")

ggplot(xplot, aes(Rsq, study_gxe)) + 
  geom_point() + 
  xlim(0,1)


xplot <- mutate_at(x, vars(matches("Frq")), as.numeric) %>% 
  pivot_longer( cols = contains("Frq"), names_to = "study_gxe", values_to = "MAF") %>% 
  mutate(MAF = 1-MAF)

ggplot(xplot, aes(MAF, study_gxe)) + 
  geom_point() + 
  xlim(0,1)




# extract chromosome 17 snp for analysis

alc_moderate_e <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_glm.rds")
alc_moderate_hit <- readRDS("~/FIGI_genotype_dosages_alcoholc_moderate.rds") %>% 
  filter(vcfid %in% alc_moderate_e$vcfid)

alc_moderate <- inner_join(alc_moderate_e, alc_moderate_hit, 'vcfid')


# try to get code for likelihood ratio and GxE meta-analysis

glm_func <- function(y) glm(outcome ~ y * alcoholc_moderate + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = alc_moderate, family = binomial)
glm_func_base <- function(y) glm(outcome ~ y + alcoholc_moderate + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = alc_moderate, family = binomial)
run_lrtest <- function(x,y) lrtest(x, y)


base_function <- glm(outcome ~ X17.27861663.G.A + alcoholc_moderate + age_ref_imp + sex + study_gxe + energytot_imp + PC1 + PC2 + PC3, data = alc_moderate, family = binomial)
int_function <- glm(outcome ~ X17.27861663.G.A * alcoholc_moderate + age_ref_imp + sex + study_gxe + energytot_imp + PC1 + PC2 + PC3, data = alc_moderate, family = binomial)
run_lrtest <- lrtest(base_function, int_function)
run_lrtest

filter(gxe, grepl("17:27861663", SNP))


base_function <- glm(alcoholc_moderate ~ X17.27861663.G.A + age_ref_imp + sex + study_gxe + energytot_imp + PC1 + PC2 + PC3, data = alc_moderate, family = binomial)
int_function <-  glm(alcoholc_moderate ~                    age_ref_imp + sex + study_gxe + energytot_imp + PC1 + PC2 + PC3, data = alc_moderate, family = binomial)
run_lrtest <- lrtest(base_function, int_function)
run_lrtest

filter(gxe, grepl("17:27861663", SNP))



# let's do meta-analysis of the GxE term... and E|G
glm_func <-      function(y) glm(alcoholc_moderate ~ y + age_ref_imp + sex + study_gxe + energytot_imp + PC1 + PC2 + PC3, data = alc_moderate, family = binomial)
glm_func_base <- function(y) glm(alcoholc_moderate ~     age_ref_imp + sex + study_gxe + energytot_imp + PC1 + PC2 + PC3, data = alc_moderate, family = binomial)
run_lrtest <- function(x,y) lrtest(x, y)







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

