#=============================================================================#
# NSAIDS - aspirin results
# 06/01/2019
# 
# Filter by Weighted Average Rsq (until another solution is devised e.g. reimputation)
# Generate Plots
# Implement 2-step methods
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())


# annotations
# (right now, I only use this df to ensure that these markers stay in plot after filtering by P < 0.05)
fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"))

# Rsq estimate based on alt allele probability
# filtered maf > 0.001, Rsq >= 0.8
rsq_filter <- readRDS("~/data/Rsq_Estimate/FIGI_RsqEstimate_chrALL.rds")

# results
gxe_all <- do.call(rbind, lapply(list.files(path = "~/data/Results/aspirin/", full.names = T, pattern = "FIGI_GxESet_aspirin_sex_age_pc3_studygxe_66485_chr"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = ":"),
         chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE) %>% 
  filter(!duplicated(ID)) # small issue with gxescan package

gxe <- gxe_all %>% 
  filter(ID %in% rsq_filter$id)


# LD Clump Results

## chiSqGxE
tmp <- do.call(rbind, lapply(list.files("~/data/Results/aspirin/clump/", full.names = T, pattern = "Plink_aspirin_ldclump_chiSqGxE_"), fread, stringsAsFactors = F))
gxe_clump_chiSqGxE <- gxe %>% 
  filter(ID %in% tmp$SNP)
rm(tmp)


## chiSqG
tmp <- do.call(rbind, lapply(list.files("~/data/Results/aspirin/clump", full.names = T, pattern = "Plink_aspirin_ldclump_chiSqG_"), fread, stringsAsFactors = F))
gxe_clump_chiSqG <- gxe %>% 
  filter(ID %in% tmp$SNP)
rm(tmp)


y <- mutate(gxe_clump_chiSqG, pval = pchisq(chiSqG, 1, lower.tail = F)) %>% 
  filter(pval < 5e-8)
x <- inner_join(gxe_clump_chiSqG, fh_annotations, by = "SNP")
x <- inner_join(y, fh_annotations, by = "SNP")


#-----------------------------------------------------------------------------#
# QQ and Manhattan Plots ----
#-----------------------------------------------------------------------------#
plot_exposure <- "aspirin"
plot_covariates <- c("age_ref_imp", "sex", "study_gxe", "PC1", "PC2", "PC3")

# Marginal G Results ----
create_qqplot(gxe_clump_chiSqG, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_DELETE")
create_manhattanplot(gxe_clump_chiSqG, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_DELETE")


# filter out fh...

x <- filter(gxe_clump_chiSqG, !SNP %in% fh_annotations$SNP)

create_manhattanplot(x, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_DELETE2")





wtf <- gwas95[[95]]
head(wtf)





load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")

table(figi$filename)
table(figi_samplefile$filename)

fucb <- filter(figi_samplefile, filename == "corect_oncoarray") %>% 
  mutate(plink_out = ifelse(outc == "", -9, 
                            ifelse(outc == "Other", -9, 
                                   ifelse(outc == "Case", 2, 1))))

table(fucb$plink_out)
out <- data.frame(IID = fucb$vcfid, 
                  FID = fucb$vcfid, 
                  pheno = fucb$plink_out)
write.table(out, file = "~/corect_oncoarray_phenotype.txt", row.names = F, quote = F, sep = '\t')


fucb <- filter(figi_samplefile, filename == "corect_oncoarray") %>% 
  mutate(sex = ifelse(sex == "Male", 1,
                      ifelse(sex == "Female", 2, 0)))



out <- data.frame(IID = fucb$vcfid, 
                  FID = fucb$vcfid, 
                  sex = fucb$sex)
write.table(out, file = "~/corect_oncoarray_sex.txt", row.names = F, quote = F, sep = '\t')




dude <- fread("~/test.ld") %>% 
  filter(SNP_B == "20:33173883:C:T")


x <- filter(gxe, grepl("7161745", SNP))

x <- filter(gxe, SNP=="14:74028117")
















###############################

# Garrett, let's see GLM for that one SNP he's talking about
# 14:74028117

chr14 <- readRDS("~/data/BinaryDosage_InfoFile/FIGI_chr14.rds")
chr14_snp <- chr14[[11]]
wtf <- dplyr::filter(chr14_snp, grepl("74028117", SNPID))

indexxx <- which(chr14_snp$SNPID == "14:74028117")


cov <- figi %>% 
  dplyr::select(vcfid, age_ref_imp, sex, studyname)
  filter(drop == 0)

  

pca <- "/home/rak/data/PCA/190729/FIGI_GwasSet_190729.eigenvec"
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")
pc30k <- fread(pca, skip = 1, col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))

cov <- figi %>%
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  filter(drop == 0,
         outc != "Other") %>% 
  dplyr::select(vcfid, outc, age_ref_imp, sex, PC1, PC2, PC3, PC4, studyname) %>% 
  filter(complete.cases(.)) %>% 
  mutate(outcome = ifelse(outc == "Control", 0, 1))

snp <- data.frame(readRDS("~/garrett_chr14.rds")) %>% 
  rownames_to_column('vcfid')

covf <- inner_join(cov, snp, 'vcfid')

ff <- glm(outcome ~ X14.74028117 + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4, data = covf, family = 'binomial')
summary(ff)

ff <- glm(outcome ~ X14.74028117 + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + studyname, data = covf, family = 'binomial')
summary(ff)
