#=============================================================================#
# 2109/09/16
# Fit main effects models with the following SNPs
# (Josh confirming MSI Positive hits on FIGI data)
#
#=============================================================================#
library(tidyverse)
library(data.table)
library(broom)
library(furrr)
library(lmtest)
rm(list = ls())

results <- fread("~/Dropbox/FIGI/Code/Misc/Josh_20190916/Top_MSI_SNPs_forFollowUpIn_corect.csv") %>% 
  mutate(SNP = paste(chr, pos, sep = ":")) %>% 
  arrange(chr, pos)

chroms <- unique(names(table(results$chr)))


# loop, output index positions as rds objects
for(chr in chroms) {
  figi <- readRDS(paste0("~/data/BinaryDosage_InfoFile/FIGI_chr", chr, ".rds"))[[11]] %>% 
    mutate(SNP = paste(Chromosome, Location, sep = ":"))
  saveRDS(which(figi$SNP %in% results$SNP), file = paste0("~/Dropbox/FIGI/Code/Misc/Josh_20190916/Josh_MSI_Hits_chr", chr, ".rds"), version = 2)
}


# output data.frame with snp information
snpinfo <- lapply(chroms, function(chr) {
  figi <- readRDS(paste0("~/data/BinaryDosage_InfoFile/FIGI_chr", chr, ".rds"))[[11]] %>% 
    filter(SNPID %in% results$SNP)
})


snpinfo_clean <- do.call(rbind, snpinfo) %>% 
  mutate(SNP = paste0("X", Chromosome, ".", Location))



#-----------------------------------------------------------------------------#
# GLM on GWAS Set
# Exclude nonEUR corect batch
#-----------------------------------------------------------------------------#
rm(list = ls())
pca <- "/home/rak/data/PCA/190729/FIGI_GwasSet_190729.eigenvec"
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")

pc <- fread(pca, skip = 1, col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))

cov <- figi %>%
  filter(drop == 0, 
         filename != "corect_oncoarray_nonEUR_reimpute") %>% 
  inner_join(pc, by = c('vcfid' = 'IID')) %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1),
         study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe))


table(cov$study_gxe, cov$outcome)
sort(unique(cov$study_gxe))
drops <- data.frame(table(cov$study_gxe, cov$outcome)) %>% 
  filter(Freq == 0)

cov <- filter(cov, !study_gxe %in% unique(drops$Var1))


# dosages
j1 <- lapply(list.files("~/Dropbox/FIGI/Code/Misc/Josh_20190916/", pattern = "out.rds", full.names = T), readRDS)
j2 <- Reduce(inner_join, j1)

cov_df <- inner_join(cov, j2, 'vcfid')

names(j2)[grepl("X", names(j2))]
head(cov_df)


results <- cov_df %>% 
  select(names(j2)[grepl("X", names(j2))]) %>% 
  purrr::map(~glm(cov_df$outcome ~ .x + cov_df$age_ref_imp+cov_df$sex+cov_df$study_gxe+cov_df$PC1+cov_df$PC2+cov_df$PC3) %>% tidy, data = cov_df, family = 'binomial') 

saveRDS(results, file = "~/Dropbox/FIGI/Code/Misc/Josh_20190916/results.rds", version = 2)

# clean results
results_wald_clean <- do.call(rbind, lapply(results, function(x) x[2,])) %>% 
  mutate(variant = names(j2)[grepl("X", names(j2))])

# write.csv(results_clean, file = "~/Dropbox/FIGI/Code/Misc/Josh_20190916/results_clean.csv", quote = F, row.names = F)



#-----------------------------------------------------------------------------#
# GLM with LRTEST
#-----------------------------------------------------------------------------#

# test with single SNP

# test_base <- glm(cov_df$outcome ~               cov_df$age_ref_imp+cov_df$sex+cov_df$study_gxe+cov_df$PC1+cov_df$PC2+cov_df$PC3, data = cov_df, family = 'binomial')
# test <-      glm(cov_df$outcome ~ X1.59209416 + cov_df$age_ref_imp+cov_df$sex+cov_df$study_gxe+cov_df$PC1+cov_df$PC2+cov_df$PC3, data = cov_df, family = 'binomial')
# 
# test_lrtest <- lrtest(test_base, test)
# test_lrtest
# 
# summary(test)


# let's apply to the markers in question
results <- cov_df %>% 
  select(names(j2)[grepl("X", names(j2))]) %>% 
  purrr::map(~lrtest(glm(cov_df$outcome ~ .x + cov_df$age_ref_imp+cov_df$sex+cov_df$study_gxe+cov_df$PC1+cov_df$PC2+cov_df$PC3),
                     glm(cov_df$outcome ~      cov_df$age_ref_imp+cov_df$sex+cov_df$study_gxe+cov_df$PC1+cov_df$PC2+cov_df$PC3)), data = cov_df, family = 'binomial') 


results_lrtest_clean <- data.frame(do.call(rbind, lapply(results, function(x) x[2,5]))) %>% 
  rownames_to_column('SNP') %>% 
  dplyr::rename(lrtest_pval = "do.call.rbind..lapply.results..function.x..x.2..5...")


# compare with wald
results_wald_clean <- fread("~/Dropbox/FIGI/Code/Misc/Josh_20190916/results_clean.csv")






#-----------------------------------------------------------------------------#
# combine all the pieces above
#-----------------------------------------------------------------------------#

results_final <- inner_join(results_wald_clean, results_lrtest_clean, by = c("variant" = "SNP")) %>% 
  inner_join(snpinfo_clean, by = c("variant" = "SNP")) %>% 
  rename(p.value.wald = p.value,
         p.value.lrtest = lrtest_pval) %>% 
  dplyr::select(variant, Chromosome, Location, Reference, Alternate, estimate, std.error, statistic, p.value.wald, p.value.lrtest)

write.csv(results_final, file = "~/Dropbox/FIGI/Code/Misc/Josh_20190916/results_clean_lrtest_pc3_124702.csv", quote = F, row.names = F)


table(cov_df$outcome)



#-----------------------------------------------------------------------------#
# allele frequency
#-----------------------------------------------------------------------------#


allelefreq <- cov_df %>% 
  select(names(j2)[grepl("X", names(j2))]) %>% 
  purrr::map(function(x) sum(x) / (2*nrow(.)))

allelefreq_clean <- data.frame(do.call(rbind, allelefreq)) %>% 
  rownames_to_column('variant') %>% 
  rename(alt_allele_freq = "do.call.rbind..allelefreq.")

write.csv(allelefreq_clean, file = "~/Dropbox/FIGI/Code/Misc/Josh_20190916/results_alt_allele_freq.csv", quote = F, row.names = F)
