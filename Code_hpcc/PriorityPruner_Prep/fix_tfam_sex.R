library(dplyr)
library(data.table)

# priority prunner has some requirements for the sex column

args <- commandArgs(trailingOnly = T)
filename <- args[1]

cov <- readRDS("FIGI_GxESet_asp_ref_sex_age_pc10_studygxe_72145_GLM.rds")

fam <- fread(paste0(filename, ".tfam")) %>%
  inner_join(cov, by = c("V1" = "vcfid")) %>%
  mutate(newsex = ifelse(sex == 0, 2, sex)) %>%
  dplyr::select(V1, V2, V3, V4, newsex, V6)
write.table(fam, file = paste0(filename, "_sexfix.tfam"),  quote = F, row.names = F, col.names = F)