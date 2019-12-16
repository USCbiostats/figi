
covariate_file <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_72269_GLM.rds")
dosages <- data.frame(readRDS("~/data/Results/aspirin/dosage/GetSNPValues_aspirin_sex_age_pc3_studygxe_index_positions_chr5_out.rds"))


# need to merge in dosage info to the EpiData
posthoc_df <- inner_join(covariate_file, dosages, by = 'vcfid') %>% 
  mutate(sex = factor(sex, levels = c(0,1), labels = c("Female", "Male")),
         aspirin = factor(aspirin, levels = c(0,1), labels = c("No", "Yes")),
         study_gxe = as.factor(study_gxe))

model1 <- glm(outcome ~ X5.40252294*aspirin + sex + age_ref_imp + study_gxe + PC1+PC2+PC3, data = posthoc_df, family = 'binomial')

model_eff <- effect("X5.40252294*aspirin",
                    model1,
                    #xlevels=xlevel_list, # might have to add more here, depending on coding of E
                    se=TRUE,
                    confidence.level=.95,
                    typical=mean)
plot(model_eff)

plot(model_eff, x.var = 'aspirin', lines = list(multilines = TRUE), xlevels = list(X5.40252294 = c(0,1,2)))


plot(predictorEffect("X5.40252294", model1, xlevels = list(X5.40252294 = c(0,1,2))), lines = list(multilines = T))
plot(predictorEffect("aspirin", model1, xlevels = list(X5.40252294 = c(0,1,2))), lines = list(multilines = TRUE))

# the other way

test <- posthoc_df %>% 
  mutate(X5.40252294 = round(X5.40252294))

model2 <- glm(outcome ~ aspirin:X5.40252294 + aspirin + X5.40252294 + sex + age_ref_imp + study_gxe + PC1+PC2+PC3, data = test, family = 'binomial')

model_eff2 <- effect("aspirin:X5.40252294",
                    model2,
                    x.var = 'aspirin',
                    xlevels=list(0,1,2), # might have to add more here, depending on coding of E
                    se=TRUE,
                    confidence.level=.95,
                    typical=mean)

plot(model_eff2, x.var = 'aspirin')

round(posthoc_df$X5.40252294)



# make sure you're doing it right 

test_aspno <- filter(test, aspirin == "No")
model2 <- glm(outcome ~ X5.40252294 + sex + age_ref_imp + study_gxe + PC1+PC2+PC3, data = test_aspno, family = 'binomial')
summary(model2)

test_aspyes <- filter(test, aspirin == "Yes")
model2 <- glm(outcome ~ X5.40252294 + sex + age_ref_imp + study_gxe + PC1+PC2+PC3, data = test_aspyes, family = 'binomial')
summary(model2)