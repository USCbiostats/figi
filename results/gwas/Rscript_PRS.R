#=============================================================================#
# Polygenic Risk Scores
# using 121 GWAS top hits (from Fred Hutch Annotation)
# excludes "Suggestive" hits
# 
# GxE with exposure variables (try with aspirin first)
# again - this is in context of aspirin GWIS (N = smaller than GWAS)
#=============================================================================#
library(tidyverse)
library(data.table)
library(purrr)
library(broom)
library(figifs)
rm(list = ls())


# let's try running glm 

j1 <- lapply(list.files("~/data/Dosages/gwas_hits_140/", pattern = "GWAS_hits_indices_chr", full.names = T), readRDS)
j2 <- Reduce(inner_join, j1) # ? why 121..?


cov_glm <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_72269_GLM.rds") %>% 
  inner_join(j2, 'vcfid')


simple <- glm(outcome ~ X8.117630683 + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = cov_glm, family = 'binomial')
summary(simple)


# run meta analysis of G main effect - take X8.117630683 as example
study_info <- readRDS("~/data/Annotations/FIGI_studygxe_info.rds")
gxe_counts <- get_counts_outcome_by_group(cov_glm, outcome, X8.117630683, study_gxe)
gxe_glm <- get_estimates_e_by_group(cov_glm, outcome, X8.117630683, group=study_gxe, age_ref_imp, sex)
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




# Create a simple version of polygenic risk score

## First, you need to get regression estimates for each of the SNPs
## 10-129 in this example

glm_g <- function(data, outcome, ...) {
  
  outcome_quo <- enquo(outcome)
  covars_quo <- enexprs(...) # it's like enquo, but it doesn't capture the environment...
  
  glmf <- function(snp) {
    glm(as.formula(paste(quo_name(outcome_quo), "~", quote(snp), "+", paste(covars_quo, collapse = " + "))), data = data, family = 'binomial')
  }

  map(data[,10:ncol(data)], ~ glmf(.x) %>% tidy, .progress = T)
}


results <- glm_g(cov_glm, outcome = 'outcome', age_ref_imp, sex, PC1, PC2, PC3, study_gxe)


## get a vector of glm estimates for the SNP

results_estimates <- do.call(cbind, lapply(results, function(x) as.numeric(x[2,2])))

cov_glm_wt <- data.frame(mapply("*", cov_glm[10:129], results_estimates))

prs <- rowSums(cov_glm_wt)


cov_prs <- cbind(cov_glm, prs)
cov_prs$prs_scaled <- scale(cov_prs$prs)

simple <- glm(outcome ~ prs + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = cov_prs, family = 'binomial')
summary(simple)

simple_scaled <- glm(outcome ~ prs_scaled + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = cov_prs, family = 'binomial')
summary(simple_scaled)


# meta analysis of the PRS result
study_info <- readRDS("~/data/Annotations/FIGI_studygxe_info.rds")
gxe_counts <- get_counts_outcome_by_group(cov_prs, outcome, prs_scaled, study_gxe)
gxe_glm <- get_estimates_e_by_group(cov_prs, outcome, prs_scaled, group=study_gxe, age_ref_imp, sex, PC1, PC2, PC3)
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


# NOW - let's see what happens when you do interaction testing between PRS and the exposure of interest

simple_gxe <- glm(outcome ~ prs_scaled*aspirin + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = cov_prs, family = 'binomial')
summary(simple_gxe) # not significant


# meta analysis of GxE 
# (also remember to add stratified tables)
# study_info <- readRDS("~/data/Annotations/FIGI_studygxe_info.rds")
# gxe_counts <- get_counts_outcome_by_group(cov_glm, outcome, X8.117630683, study_gxe)
# gxe_glm <- get_estimates_gxe_by_group(cov_glm, outcome, aspirin, study_gxe, X8.117630683, age_ref_imp, sex)
# gxe_meta <- dplyr::bind_cols(gxe_counts, gxe_glm) %>%
#   inner_join(., study_info, by = "study_gxe")
# 
# results_meta <- meta::metagen(estimate,
#                               std.error,
#                               data=gxe_meta,
#                               studlab=paste(study_gxe),
#                               comb.fixed = FALSE,
#                               comb.random = TRUE,
#                               method.tau = "SJ",
#                               hakn = TRUE,
#                               prediction=TRUE,
#                               sm="OR",
#                               byvar = study_design)
# 
# meta::forest(results_meta,
#              layout = "JAMA",
#              text.predict = "95% CI",
#              col.predict = "black",
#              leftcols = c("studlab", "Control", "Case", "N", "effect", "ci", "w.random"),
#              digits.addcols=0)




