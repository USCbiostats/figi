#=============================================================================#
# FIGI Analysis 02/26/2019
# FIGI Analysis 03/13/2019
# FIGI Analysis 05/16/2019
# FIGI Analysis 08/08/2019 - reimputed UKB + reach
# FIGI Analysis 12/05/2019
#
# simple R object of GxE dataset
# without and without non-EUR
#
# sometimes sex/outcome need to be numeric, sometimes factors.. 
# in figi_gxe, outcome and sex are numeric 
#=============================================================================#
library(tidyverse)
library(data.table)
rm(list = ls())
figi_gwas <- readRDS("~/data/FIGI_EpiData_rdata/FIGI_HRC_v2.3_GWAS.rds")


figi_gxe <- figi_gwas %>% 
  dplyr::filter(gxe == 1) %>%
  dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
                sex = ifelse(sex == "Female", 0, 1))

figi_gxe_EUR <- figi_gwas %>% 
  dplyr::filter(gxe == 1,
                EUR_subset == 1) %>%
  dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
                sex = ifelse(sex == "Female", 0, 1))


saveRDS(figi_gxe,     file = '~/data/GxEScanR_PhenotypeFiles/FIGI_HRC_v2.3_GXE.rds', version = 2)
saveRDS(figi_gxe_EUR, file = '~/data/GxEScanR_PhenotypeFiles/FIGI_HRC_v2.3_GXE_EUR.rds', version = 2)




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




