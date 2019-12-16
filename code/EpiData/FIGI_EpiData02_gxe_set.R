#=============================================================================#
# FIGI Analysis 02/26/2019
# FIGI Analysis 03/13/2019
# FIGI Analysis 05/16/2019
# FIGI Analysis 08/08/2019 - reimputed UKB + reach
#
# simple R object of GxE dataset
# without and without non-EUR
#=============================================================================#
library(tidyverse)
library(data.table)
rm(list = ls())
pca <- "/home/rak/data/PCA/190729/FIGI_GxESet_190729.eigenvec"
pc <- fread(pca, skip = 1, col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")


#-----------------------------------------------------------------------------#
# GxE Set asp_ref
#-----------------------------------------------------------------------------#
cov_gxe <- figi %>%
  filter(drop == 0 & gxe == 1) %>% 
  inner_join(pc, by = c('vcfid' = 'IID'))

table(cov_gxe$EUR_subset)

cov_gxe_EUR <- cov_gxe %>% 
  filter(EUR_subset == 1)

table(cov_gxe_EUR$race_self)

saveRDS(cov_gxe,     file = '~/data/GxEScanR_PhenotypeFiles/FIGI_gxeset.rds', version = 2)
saveRDS(cov_gxe_EUR, file = '~/data/GxEScanR_PhenotypeFiles/FIGI_gxeset_EUR.rds', version = 2)




#-----------------------------------------------------------------------------#
# TEMPORARY - NO PCs
#-----------------------------------------------------------------------------#
library(tidyverse)
library(data.table)
rm(list = ls())
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_191205.RData")

cov_gxe <- figi %>%
  filter(drop == 0 & gxe == 1)

table(cov_gxe$EUR_subset)

cov_gxe_EUR <- cov_gxe %>%
  filter(EUR_subset == 1)

table(cov_gxe_EUR$race_self)

saveRDS(cov_gxe,     file = '~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset.rds', version = 2)
saveRDS(cov_gxe_EUR, file = '~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_EUR.rds', version = 2)




