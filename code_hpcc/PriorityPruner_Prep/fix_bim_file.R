library(tidyverse)
library(data.table)

# need this step to ensure that every batch has same SNP IDs

args <- commandArgs(trailingOnly = T)
filename <- args[1]

x <- fread(paste0(filename, ".bim")) %>%
    mutate(V2 = paste0(V1, ":", V4, "_", V6, "_", V5))
write.table(x, file = paste0(filename, ".bim"), quote = F, row.names = F, col.names = F, sep = "\t")