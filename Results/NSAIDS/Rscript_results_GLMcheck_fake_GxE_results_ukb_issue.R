#=======================================================================================#
# Perform check using GLM
# of the 10 'fake' GxE GWIS hits 
# (that turned out to be discrepant allele frequencies in UKB)
#=======================================================================================#
library(lmtest)
library(tidyverse)
library(data.table)

# ............... after extracting from bdose files ..............#
# (see the results_asp_ref_qq_manhattan.R script)

# covariate file
cov <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc10_studygxe_72145_GLM.rds")

for(chr in c(2,3,6,7,10,16,17,18,21,22)) {
  x <- as.data.frame(readRDS(paste0("./qc/GxEScanR_GxE_sig_loci_extract_chr", chr, ".rds"))) %>% 
    rownames_to_column(var = 'vcfid')
  cov <- inner_join(cov, x, by = 'vcfid')
}



# Let's replicate the results from GxEScanR
# Remember it's lrt not wald
glm_func <- function(y) glm(outcome ~ y * asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + study_gxe, data = cov, family = binomial)
glm_func_base <- function(y) glm(outcome ~ y + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + study_gxe, data = cov, family = binomial)

basemodel <- map(cov[,17:26], function(x) glm_func_base(x))
intmodel <-           map(cov[,17:26], function(x) glm_func(x))
run_lrtest <- function(x,y) lrtest(x, y)
results_gxe <- mapply(run_lrtest, basemodel, intmodel, SIMPLIFY = F)

results_gxe_clean <- do.call(rbind, results_gxe) %>% 
  tibble::rownames_to_column('SNP') %>% 
  filter(!is.na(Chisq)) %>% 
  mutate(SNP = gsub('.{2}$', '', SNP))

z <- inner_join(results_GxE_filter, results_gxe_clean, by = 'SNP') %>% 
  mutate(betaGxE_gxescan = betaGxE, 
         P.GxE_gxescan = P, 
         P.GxE_glm = `Pr(>Chisq)`) %>% 
  dplyr::select(SNP, Subjects, Cases, betaGxE_gxescan, P.GxE_gxescan, P.GxE_glm)

plot(-log10(z$P.GxE_gxescan), -log10(z$P.GxE_glm))
abline(a = 0, b = 1) #OK. 


# what if you DIDN'T adjust for study... 
# (P values become more extreme)
glm_func <- function(y) glm(outcome ~ y * asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = cov, family = binomial)
glm_func_base <- function(y) glm(outcome ~ y + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = cov, family = binomial)

basemodel <- map(cov[,17:26], function(x) glm_func_base(x))
intmodel <-           map(cov[,17:26], function(x) glm_func(x))
run_lrtest <- function(x,y) lrtest(x, y)
results_gxe <- mapply(run_lrtest, basemodel, intmodel, SIMPLIFY = F)

results_gxe_clean <- do.call(rbind, results_gxe) %>% 
  tibble::rownames_to_column('SNP') %>% 
  filter(!is.na(Chisq)) %>% 
  mutate(SNP = gsub('.{2}$', '', SNP))

z <- inner_join(results_GxE_filter, results_gxe_clean, by = 'SNP') %>% 
  mutate(betaGxE_gxescan = betaGxE, 
         P.GxE_gxescan = P, 
         P.GxE_glm = `Pr(>Chisq)`) %>% 
  dplyr::select(SNP, Subjects, Cases, betaGxE_gxescan, P.GxE_gxescan, P.GxE_glm)

plot(-log10(z$P.GxE_gxescan), -log10(z$P.GxE_glm))
abline(a = 0, b = 1)
