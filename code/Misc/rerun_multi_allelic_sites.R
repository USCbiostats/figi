#=============================================================================#
# Andre Kim
# 01/10/2020
#
# problem - multi-allelic sites were being duplicated in gxescan results
# such that we'd potentially miss results from the third alternate allele marker
#
#
# output index positions of multi-allelic sites, and run gxescan for those
# sites only, then merge these results with the results from the original
# scans. This way, you'll have a complete set without having to re-run 
# the whole thing
#=============================================================================#
library(tidyverse)
library(data.table)

x <- fread("~/data/HRC/HRC.r1-1.GRCh37.chr22.txt") %>% 
  mutate(id = paste0(V1, ":", V2, ":", V3, ":", V4))
hrc <- readRDS("~/data/BinaryDosage_InfoFile/FIGI_snpid_fix_chr22.rds")$SNPs

dups <- x[which(duplicated(x$V2)), ]
x_dups <- filter(x, V2 %in% dups$V2)

ind <- which(hrc$SNPID %in% x_dups$id)
saveRDS()


# write and run function



wrap <- function(chr) {
  
  x <- fread(paste0("~/data/HRC/HRC.r1-1.GRCh37.chr", chr, ".txt")) %>% 
    mutate(id = paste0(V1, ":", V2, ":", V3, ":", V4))
  hrc <- readRDS(paste0("~/data/BinaryDosage_InfoFile/FIGI_snpid_fix_chr", chr, ".rds"))$SNPs
  
  dups <- x[which(duplicated(x$V2)), ]
  x_dups <- filter(x, V2 %in% dups$V2)
  
  ind <- which(hrc$SNPID %in% x_dups$id)
  saveRDS(ind, file = paste0("/home/rak/tmp_multiallelics/multiallelic_locations/hrc_multiallelic_sites_index_chr", chr, ".rds"), version = 2)
  
}

for(x in 1:21) {wrap(x)}
