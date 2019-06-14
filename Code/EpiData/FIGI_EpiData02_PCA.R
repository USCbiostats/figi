#-----------------------------------------------------------
# FIGI Analysis 10/26/2018
#
# Get GWAS and GxE sets for PCA re-calculation
# output FID \t IID (for plink subsetting)
# (vcfid \t vcfid)
#-----------------------------------------------------------
library(tidyverse)
library(data.table)
rm(list = ls())

load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190424.RData")


# ------ GWAS Set ------
# N = 138,572
gwas_set <- filter(figi, drop == 0) %>% 
  mutate(FID = vcfid,
         IID = vcfid) %>% 
  dplyr::select(FID, IID)

any(duplicated(gwas_set$vcfid))
write.table(gwas_set, file = "~/FIGI_PCA_GwasSet_190430.txt", quote = F, row.names = F, col.names = F, sep = '\t')



# ------ GxE Set ------
# N = 102,792
gxe_set <- filter(figi, drop == 0, gxe == 1) %>% 
  mutate(FID = vcfid,
         IID = vcfid) %>% 
  dplyr::select(FID, IID)

any(duplicated(gxe_set$vcfid))
write.table(gxe_set, file = "~/FIGI_PCA_GxESet_190430.txt", quote = F, row.names = F, col.names = F, sep = '\t')








#-----------------------------------------------------------------------------#



# ------ GWAS Set noUKB (tmp) ------
# N = 138,572
gwas_set <- filter(figi, drop == 0, study_gxe != "UKB_1") %>% 
  mutate(FID = vcfid,
         IID = vcfid) %>% 
  dplyr::select(FID, IID)

any(duplicated(gwas_set$vcfid))
write.table(gwas_set, file = "~/FIGI_PCA_GwasSet_noUKB_190430.txt", quote = F, row.names = F, col.names = F, sep = '\t')



# ------ GxE Set noUKB (tmp) ------
# N = 102,792
gxe_set <- filter(figi, drop == 0, gxe == 1, study_gxe != "UKB_1") %>% 
  mutate(FID = vcfid,
         IID = vcfid) %>% 
  dplyr::select(FID, IID)

any(duplicated(gxe_set$vcfid))
write.table(gxe_set, file = "~/FIGI_PCA_GxESet_noUKB_190430.txt", quote = F, row.names = F, col.names = F, sep = '\t')




