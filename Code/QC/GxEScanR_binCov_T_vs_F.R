#=============================================================================#
# GxEScanR
# 
# compare options binCov=T vs binCov=F
# (difference should be small, but good to check statistics)
#=============================================================================#
library(data.table)
library(tidyverse)

figi_chr22_binCovT <- fread("~/data/Results/NSAIDS/asp_ref_190518/results_GxE_asp_ref_sex_age_pc3_studygxe_72820_binCovT_chr22.out")
figi_chr22_binCovF <- fread("~/results_GxE_asp_ref_sex_age_pc3_studygxe_72820_binCovF_chr22.out")

plot(figi_chr22_binCovT$chiSqGE, figi_chr22_binCovF$chiSqGE)
