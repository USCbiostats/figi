#=============================================================================#
# does changing reference group change estimates or p values?
# NO 
#=============================================================================#
figi_gwas <- readRDS("~/data/FIGI_EpiData_rdata/FIGI_HRC_v2.3_GWAS.rds")
figi_gxe <- dplyr::filter(figi_gwas, gxe == 1, EUR_subset == 1)


g <- readRDS("~/huygue_gwas_140_tophitonly_10_101315166_TMP1.rds")

work <- inner_join(figi_gxe, g, 'vcfid') %>% 
  filter(outcome != "Other") %>% 
  mutate(e1 = factor(asp_ref, levels = c("No","Yes")), 
         e2 = factor(asp_ref, levels = c("Yes","No")), 
         outcome = fct_drop(outcome)) 

table(work$e1)
levels(work$e1)
table(work$e2)
levels(work$e2)
levels(work$outcome)



model1 <- glm(outcome ~ e1 + X10.101315166 + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = work, family = 'binomial')
summary(model1)
#e1Yes                    2.752e-01  1.717e-02  16.027  < 2e-16 ***

model2 <- glm(outcome ~ e2 + X10.101315166 + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = work, family = 'binomial')
summary(model2)
#e2No                    -2.752e-01  1.717e-02 -16.027  < 2e-16 *** INVERSE




# how about for gxe
model1 <- glm(outcome ~ e1 * X10.101315166 + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = work, family = 'binomial')
summary(model1)
#e1Yes                    2.589e-01  2.190e-02  11.825  < 2e-16 ***
#e1Yes:X10.101315166      3.524e-02  2.955e-02   1.193 0.232957

model2 <- glm(outcome ~ e2 * X10.101315166 + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = work, family = 'binomial')
summary(model2)
#e2No:X10.101315166      -3.524e-02  2.955e-02  -1.193 0.232957 INVERSE
#e2No                    -2.589e-01  2.190e-02 -11.825  < 2e-16 ***




# make sure of this by running lrtest
library(lmtest)

base1 <- glm(outcome ~ e1 + X10.101315166 + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = work, family = 'binomial')
model1 <- glm(outcome ~ e1 * X10.101315166 + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = work, family = 'binomial')
lrtest(base1, model1)

# Likelihood ratio test
# 
# Model 1: outcome ~ e1 + X10.101315166 + age_ref_imp + sex + study_gxe + 
#   PC1 + PC2 + PC3
# Model 2: outcome ~ e1 * X10.101315166 + age_ref_imp + sex + study_gxe + 
#   PC1 + PC2 + PC3
# #Df LogLik Df  Chisq Pr(>Chisq)
# 1  64 -45122                     
# 2  65 -45121  1 1.4233     0.2329

base2 <- glm(outcome ~ e2 + X10.101315166 + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = work, family = 'binomial')
model2 <- glm(outcome ~ e2 * X10.101315166 + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = work, family = 'binomial')
lrtest(base2, model2)


# Likelihood ratio test
# 
# Model 1: outcome ~ e2 + X10.101315166 + age_ref_imp + sex + study_gxe + 
#   PC1 + PC2 + PC3
# Model 2: outcome ~ e2 * X10.101315166 + age_ref_imp + sex + study_gxe + 
#   PC1 + PC2 + PC3
# #Df LogLik Df  Chisq Pr(>Chisq)
# 1  64 -45122                     
# 2  65 -45121  1 1.4233     0.2329








# ---- what about G|E type stuff

model3 <- glm(e1 ~ X10.101315166 + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = work, family = 'binomial')
summary(model3)
#X10.101315166           -0.0084113  0.0140717  -0.598 0.550013    

model4 <- glm(e2 ~ X10.101315166 + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = work, family = 'binomial')
summary(model4)
#X10.101315166            8.411e-03  1.407e-02   0.598 0.550013    