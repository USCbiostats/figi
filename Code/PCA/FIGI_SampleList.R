#=============================================================================#
#
# get sample lists for calculating PC/IBD
#
#=============================================================================#
library(tidyverse)
library(data.table)
rm(list = ls())
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")


#-----------------------------------------------------------------------------#
# PCA GwasSet
#-----------------------------------------------------------------------------#
cov <- figi %>%
	dplyr::filter(drop == 0)

write.table(cov[, c('vcfid', 'vcfid')], file = "~/data/PCA/190729/FIGI_PCA_GwasSet_190729.txt", quote = F, row.names = F, col.names = F, sep = '\t')

#-----------------------------------------------------------------------------#
# PCA GxESet
#-----------------------------------------------------------------------------#
cov <- figi %>%
	dplyr::filter(drop == 0 & gxe == 1)

write.table(cov[, c('vcfid', 'vcfid')], file = "~/data/PCA/190729/FIGI_PCA_GxESet_190729.txt", quote = F, row.names = F, col.names = F, sep = '\t')