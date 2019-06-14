#=============================================================================#
# Investigate source of induced GxE in the context of 
# flipped allelic dosages in UKB data vs all other samples
#=============================================================================#

# Covariate file
cov <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc10_studygxe_72145_GLM.rds")

# Dosages

# start with chr3
dos <- readRDS("./gxe_problems/GxEScanR_GxE_sig_loci_extract_chr3.rds") 
dos <- data.frame(dos) %>% 
  rownames_to_column("vcfid") %>% 
  rename(dose = "X3.98981063")
names(dos)

# do the deed
df <- inner_join(cov, dos, by = "vcfid")
names(df)


# GLM (GxE) overall

gxe_all <- glm(outcome ~ dose*asp_ref+age_ref_imp+sex+study_gxe+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df, family = 'binomial')
coef(summary(gxe_all))

# don't adjust for study when doing UKB only
gxe_ukb <- glm(outcome ~ dose*asp_ref+age_ref_imp+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df[which(df$study_gxe=="UKB_1"),], family = 'binomial')
coef(summary(gxe_ukb)) # dose:asp_ref p = 3.51e-2

gxe_nonukb <- glm(outcome ~ dose*asp_ref+age_ref_imp+sex+study_gxe+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df[which(df$study_gxe!="UKB_1"),], family = 'binomial')
coef(summary(gxe_nonukb)) # dose:asp_ref p = 5.31e-1

# what if you create variable UKB vs nonUKB
df <- df %>% 
  mutate(ukb = ifelse(study_gxe == "UKB_1", 1, 0))
names(df)


gxe_all_mod <- glm(outcome ~ dose*asp_ref+age_ref_imp+sex+study_gxe+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ukb, data = df, family = 'binomial')
coef(summary(gxe_all_mod))
