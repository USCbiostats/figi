#=============================================================================#
# FIGI Analysis 02/25/2019
#
# Compare UKB dosage vs those of others (ukbiobank vs rest indicator variable, adj for PCs)
#
#=============================================================================#
library(tidyverse)
library(data.table)

rm(list = ls())
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_181120.RData")

#------ GxE Set - asp_ref (N = 102,792 --> 72,438) ------
pc30k <- fread("~/data/PrincipalComponents/FIGI_GwasSet.eigenvec", skip = 1, 
               col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))

cov <- Epi %>%
  filter(drop == 0) %>% 
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  mutate(ukb = ifelse(studyname == "UKB_1", 1, 0)) %>% 
  dplyr::select(vcfid, ukb, paste0(rep("PC", 10), seq(1,10)))


table(cov$ukb)
head(cov)

# save cov file as rds + text
saveRDS(cov, file = '~/FIGI_GWASSet_ukb_allelefreq_check.rds')


