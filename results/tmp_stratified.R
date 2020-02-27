


rm(list = ls())
covariate_file <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_redmeatqc2_basic_covars_glm.rds")
dosages <- readRDS("~/FIGI_genotype_dosages_stern_gxe.rds")


# need to merge in dosage info to the EpiData
posthoc_df <- inner_join(covariate_file, dosages, by = 'vcfid') %>% 
  

model1 <- glm(outcome ~ X18.46451873.G.A*redmeatqc2 + sex + age_ref_imp + energytot_imp + study_gxe + PC1+PC2+PC3, data = posthoc_df, family = 'binomial')

model_eff <- effect("X18.46451873.G.A*redmeatqc2",
                    model1,
                    #xlevels=xlevel_list, # might have to add more here, depending on coding of E
                    se=TRUE,
                    confidence.level=.95,
                    typical=mean)
plot(model_eff)

plot(model_eff, x.var = 'redmeatqc2', lines = list(multilines = TRUE), xlevels = list(X18.46451873.G.A = c(0,1,2)))


plot(predictorEffect("X18.46451873.G.A", model1, xlevels = list(X18.46451873.G.A = c(0,1,2))), lines = list(multilines = T))
plot(predictorEffect("redmeatqc2", model1, xlevels = list(X18.46451873.G.A = c(0,1,2))), lines = list(multilines = TRUE))

# the other way

test <- posthoc_df %>% 
  mutate(X18.46451873.G.A = round(X18.46451873.G.A))

model2 <- glm(outcome ~ redmeatqc2:X18.46451873.G.A + redmeatqc2 + X18.46451873.G.A + sex + age_ref_imp + energytot_imp + study_gxe + PC1+PC2+PC3, data = test, family = 'binomial')

model_eff2 <- effect("redmeatqc2:X18.46451873.G.A",
                     model2,
                     x.var = 'redmeatqc2',
                     xlevels=list(0,1,2,3), # might have to add more here, depending on coding of E
                     se=TRUE,
                     confidence.level=.95,
                     typical=mean)

plot(model_eff2, x.var = 'redmeatqc2')

round(posthoc_df$X5.40252294)



# make sure you're doing it right 

test_aspno <- filter(test, aspirin == "No")
model2 <- glm(outcome ~ X5.40252294 + sex + age_ref_imp + study_gxe + PC1+PC2+PC3, data = test_aspno, family = 'binomial')
summary(model2)

test_aspyes <- filter(test, aspirin == "Yes")
model2 <- glm(outcome ~ X5.40252294 + sex + age_ref_imp + study_gxe + PC1+PC2+PC3, data = test_aspyes, family = 'binomial')
summary(model2)






# can't get it working nice
test <- posthoc_df %>% 
  mutate(X18.46451873.G.A = round(X18.46451873.G.A),
         chr18_46451873 = as.factor(X18.46451873.G.A),
         redmeatqc2 = as.factor(redmeatqc2))
table(test$chr18_46451873)

model2 <- glm(outcome ~ redmeatqc2:chr18_46451873 + redmeatqc2 + chr18_46451873 + sex + age_ref_imp + energytot_imp + study_gxe + PC1+PC2+PC3, data = test, family = 'binomial')

model_eff2 <- effect("redmeatqc2:chr18_46451873",
                     model2,
                     x.var = 'redmeatqc2',
                     xlevels=list(0,1,2,3), # might have to add more here, depending on coding of E
                     se=TRUE,
                     confidence.level=.95,
                     typical=mean)

png(paste0("~/Dropbox/redmeatqc2_categorical_stratified_plot_chr18_46451873.png"), height = 400, width = 600)
plot(model_eff2, x.var = 'redmeatqc2')
dev.off()





# can't get it working nice
test <- posthoc_df %>% 
  mutate(X18.46451873.G.A = round(X18.46451873.G.A),
         chr18_46451873 = as.factor(X18.46451873.G.A))
table(test$chr18_46451873)

model2 <- glm(outcome ~ redmeatqc2:chr18_46451873 + redmeatqc2 + chr18_46451873 + sex + age_ref_imp + energytot_imp + study_gxe + PC1+PC2+PC3, data = test, family = 'binomial')

model_eff2 <- effect("redmeatqc2:chr18_46451873",
                     model2,
                     x.var = 'redmeatqc2',
                     xlevels=list(0,1,2,3), # might have to add more here, depending on coding of E
                     se=TRUE,
                     confidence.level=.95,
                     typical=mean)


png(paste0("~/Dropbox/redmeatqc2_continuous_stratified_plot_chr18_46451873.png"), height = 400, width = 600)
plot(model_eff2, x.var = 'redmeatqc2')
dev.off()





# can't get it working nice
test <- posthoc_df %>% 
  mutate(X18.46451873.G.A = round(X18.46451873.G.A),
         chr18_46451873 = as.factor(X18.46451873.G.A))
table(test$chr18_46451873)

model2 <- glm(outcome ~ redmeatqc2:X18.46451873.G.A + redmeatqc2 + X18.46451873.G.A + sex + age_ref_imp + energytot_imp + study_gxe + PC1+PC2+PC3, data = posthoc_df, family = 'binomial')

model_eff2 <- effect("redmeatqc2:X18.46451873.G.A",
                     model2,
                     x.var = 'redmeatqc2',
                     xlevels=list(X18.46451873.G.A = c(0,1,2)), # might have to add more here, depending on coding of E
                     se=TRUE,
                     confidence.level=.95,
                     typical=mean)


png(paste0("~/Dropbox/redmeatqc2_continuous_stratified_plot_chr18_46451873.png"), height = 400, width = 600)
plot(model_eff2, x.var = 'redmeatqc2')
dev.off()






































# can't get it working nice
test <- posthoc_df %>% 
  mutate(chr18_46451873 = X18.46451873.G.A,
         redmeatqc2 = as.factor(redmeatqc2))

model2 <- glm(outcome ~ redmeatqc2:chr18_46451873 + redmeatqc2 + chr18_46451873 + sex + age_ref_imp + energytot_imp + study_gxe + PC1+PC2+PC3, data = test, family = 'binomial')

model_eff2 <- effect("redmeatqc2:chr18_46451873",
                     model2,
                     x.var = 'redmeatqc2',
                     xlevels=list(chr18_46451873 = c(0,1,2)), # might have to add more here, depending on coding of E
                     se=TRUE,
                     confidence.level=.95,
                     typical=mean)

png(paste0("~/Dropbox/redmeatqc2_categorical_stratified_plot_chr18_46451873.png"), height = 400, width = 600)
plot(model_eff2, x.var = 'redmeatqc2')
dev.off()





# can't get it working nice
test <- posthoc_df %>% 
  mutate(chr18_46451873 = X18.46451873.G.A)

model2 <- glm(outcome ~ redmeatqc2:chr18_46451873 + redmeatqc2 + chr18_46451873 + sex + age_ref_imp + energytot_imp + study_gxe + PC1+PC2+PC3, data = test, family = 'binomial')

model_eff2 <- effect("redmeatqc2:chr18_46451873",
                     model2,
                     x.var = 'redmeatqc2',
                     xlevels=list(chr18_46451873 = c(0,1,2)), # might have to add more here, depending on coding of E
                     se=TRUE,
                     confidence.level=.95,
                     typical=mean)


png(paste0("~/Dropbox/redmeatqc2_continuous_stratified_plot_chr18_46451873.png"), height = 400, width = 600)
plot(model_eff2, x.var = 'redmeatqc2')
dev.off()







