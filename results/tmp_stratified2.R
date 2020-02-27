


rm(list = ls())
covariate_file <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_fruitqc2_basic_covars_glm.rds")
dosages <- readRDS("~/FIGI_genotype_dosages_stern_gxe.rds")


# need to merge in dosage info to the EpiData
posthoc_df <- inner_join(covariate_file, dosages, by = 'vcfid')
  

# can't get it working nice
test <- posthoc_df %>% 
  mutate(chr14_74029409 = X14.74029409.C.T,
         fruitqc2 = as.factor(fruitqc2))

model2 <- glm(outcome ~ fruitqc2:chr14_74029409 + fruitqc2 + chr14_74029409 + sex + age_ref_imp + energytot_imp + study_gxe + PC1+PC2+PC3, data = test, family = 'binomial')

model_eff2 <- effect("fruitqc2:chr14_74029409",
                     model2,
                     x.var = 'fruitqc2',
                     xlevels=list(chr14_74029409 = c(0,1,2)), # might have to add more here, depending on coding of E
                     se=TRUE,
                     confidence.level=.95,
                     typical=mean,
                     latent = T)

png(paste0("~/Dropbox/fruitqc2_categorical_stratified_plot_chr14_74029409.png"), height = 400, width = 600)
plot(model_eff2, x.var = 'fruitqc2')
dev.off()





# can't get it working nice
test <- posthoc_df %>% 
  mutate(chr14_74029409 = X14.74029409.C.T)

model2 <- glm(outcome ~ fruitqc2:chr14_74029409 + fruitqc2 + chr14_74029409 + sex + age_ref_imp + energytot_imp + study_gxe + PC1+PC2+PC3, data = test, family = 'binomial')

model_eff2 <- effect("fruitqc2:chr14_74029409",
                     model2,
                     x.var = 'fruitqc2',
                     xlevels=list(chr14_74029409 = c(0,1,2)), # might have to add more here, depending on coding of E
                     se=TRUE,
                     confidence.level=.95,
                     typical=mean,
                     latent = T)


png(paste0("~/Dropbox/fruitqc2_continuous_stratified_plot_chr14_74029409.png"), height = 400, width = 600)
plot(model_eff2, x.var = 'fruitqc2')
dev.off()












# can't get it working nice
test <- posthoc_df %>% 
  mutate(chr14_74029049 = X14.74029049.G.C,
         fruitqc2 = as.factor(fruitqc2))

model2 <- glm(outcome ~ fruitqc2:chr14_74029049 + fruitqc2 + chr14_74029049 + sex + age_ref_imp + energytot_imp + study_gxe + PC1+PC2+PC3, data = test, family = 'binomial')

model_eff2 <- effect("fruitqc2:chr14_74029049",
                     model2,
                     x.var = 'fruitqc2',
                     xlevels=list(chr14_74029049 = c(0,1,2)), # might have to add more here, depending on coding of E
                     se=TRUE,
                     confidence.level=.95,
                     typical=mean)

png(paste0("~/Dropbox/fruitqc2_categorical_stratified_plot_chr14_74029049.png"), height = 400, width = 600)
plot(model_eff2, x.var = 'fruitqc2')
dev.off()





# can't get it working nice
test <- posthoc_df %>% 
  mutate(chr14_74029049 = X14.74029049.G.C)

model2 <- glm(outcome ~ fruitqc2:chr14_74029049 + fruitqc2 + chr14_74029049 + sex + age_ref_imp + energytot_imp + study_gxe + PC1+PC2+PC3, data = test, family = 'binomial')

model_eff2 <- effect("fruitqc2:chr14_74029049",
                     model2,
                     x.var = 'fruitqc2',
                     xlevels=list(chr14_74029049 = c(0,1,2)), # might have to add more here, depending on coding of E
                     se=TRUE,
                     confidence.level=.95,
                     typical=mean)


png(paste0("~/Dropbox/fruitqc2_continuous_stratified_plot_chr14_74029049.png"), height = 400, width = 600)
plot(model_eff2, x.var = 'fruitqc2')
dev.off()


