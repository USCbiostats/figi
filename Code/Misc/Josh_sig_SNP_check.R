#=============================================================================#
# Testing GxG interaction for Josh Millstein
# rs11165845   1 97819405
# rs7777031   7  8149596
#=============================================================================#


chr1 <- readRDS("/home/rak/data/FIGI_BDose_IndexFiles/FIGI_chr1.rds")
chr1_index <- which(chr1[[11]]$SNPID == "1:97819405")


chr7 <- readRDS("/home/rak/data/FIGI_BDose_IndexFiles/FIGI_chr7.rds")
chr7_index <- which(chr7[[11]]$SNPID == "7:8149596")




#-----------------------------------------------------------------------------#
# GWAS Set
# Exclude nonEUR corect batch.. 
#-----------------------------------------------------------------------------#
rm(list = ls())
pca <- "/home/rak/data/PCA/190729/FIGI_GwasSet_190729.eigenvec"
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")

pc30k <- fread(pca, skip = 1, col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))

cov <- figi %>%
  filter(drop == 0, 
         filename != "corect_oncoarray_nonEUR_reimpute") %>% 
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1),
         study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe))
  

table(cov$study_gxe, cov$outcome)
sort(unique(cov$study_gxe))
drops <- data.frame(table(cov$study_gxe, cov$outcome)) %>% 
  filter(Freq == 0)

cov <- filter(cov, !study_gxe %in% unique(drops$Var1))


# incorporate genotypes of interest
chr1 <- data.frame(readRDS("~/chr1_97819405.rds")) %>% 
  rownames_to_column('vcfid')

chr7 <- data.frame(readRDS("~/chr7_8149596.rds")) %>% 
  rownames_to_column('vcfid')

covf <- cov %>% 
  inner_join(chr1, 'vcfid') %>% 
  inner_join(chr7, 'vcfid')


# fit models (single variant)
model1 <- glm(outcome ~ X1.97819405 + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = covf, family = 'binomial')
summary(model1)

model2 <- glm(outcome ~ X7.8149596 + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = covf, family = 'binomial')
summary(model2)

model3 <- glm(outcome ~ X1.97819405*X7.8149596 + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = covf, family = 'binomial')
summary(model3)






model1 <- glm(outcome ~ X1.97819405 + age_ref_imp + sex + study_gxe , data = covf, family = 'binomial')
summary(model1)

model2 <- glm(outcome ~ X7.8149596 + age_ref_imp + sex + study_gxe , data = covf, family = 'binomial')
summary(model2)

model3 <- glm(outcome ~ X1.97819405*X7.8149596 + age_ref_imp + sex + study_gxe , data = covf, family = 'binomial')
summary(model3)





# create hardcalls, dominant coding 
# start with <=0.9, >0.9

covfh <- covf %>% 
  mutate(X1.97819405h = ifelse(X1.97819405 <= 0.9, 0, 1),
         X7.8149596h = ifelse(X7.8149596 <= 0.9, 0, 1))


model1 <- glm(outcome ~ X1.97819405h + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = covfh, family = 'binomial')
summary(model1)

model2 <- glm(outcome ~ X7.8149596h + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = covfh, family = 'binomial')
summary(model2)

model3 <- glm(outcome ~ X1.97819405h*X7.8149596h + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = covfh, family = 'binomial')
summary(model3)


