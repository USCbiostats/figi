#=============================================================================#
# Extract SNPs
#
# in trying to make this a general function, maybe the input
# should always be a vector of SNPs in the form
# chr:bp:ref:alt
#=============================================================================#

snps <- fread("~/Dropbox/FIGI/Results/aspirin/aspirin_bin2_hit_ld_snplist.txt", header = F, col.names = "ID")

# using snps, get index position
# in this case, chromosome 5 only

chr <- 5
figi <- readRDS(paste0("~/data/BinaryDosage_InfoFile/FIGI_chr", chr, ".rds"))[[11]] %>% 
  mutate(ID = paste0(Chromosome, ":", Location, ":", Reference, ":", Alternate))

saveRDS(which(figi$ID %in% snps$ID), file = "~/Dropbox/FIGI/Results/aspirin/aspirin_bin2_hit_ld_snplist_index.rds", version = 2)

