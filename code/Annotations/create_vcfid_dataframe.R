#=============================================================================#
# Create data.frame containing vcfids of EUR CONTROLS
#=============================================================================#

library(tidyverse)
library(data.table)

load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_191205.RData")

table(figi$outc)
x <- filter(figi, drop == 0, EUR_subset == 1, outc == "Control") %>% 
  dplyr::select(vcfid)

saveRDS(x, file = "/home/rak/data/Annotations/FIGI_controls_vcfid_73598.rds")
