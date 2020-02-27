#=============================================================================#
# FIGI - Extract Dosages from BinaryDosage files
# 02/04/2020
#
# Since I'm using a different info file with unique SNP names, maybe
# I can extract dosages using chr:bp:ref:alt names directly 
# to simplify one-off type tasks
#
# keep number of SNPs to reasonable level (< 20) 
#=============================================================================#
library(BinaryDosage)

# arguments
args <- commandArgs(trailingOnly=T)
filename <- args[1]
snps <- args[-1]

# process SNPs into a list of vectors, separate by chromosome
tmp1 <- sapply(snps, function(x) strsplit(x, split = ":"))
tmp2 <- as.numeric(sapply(tmp1, function(x) x[1]))
snps_list <- split(snps, tmp2)

GetSNPValues_wrapper <- function(vector_of_snps) {
  # prepare
  chr <- strsplit(vector_of_snps[1], split = ":")[[1]][1] # get chromosome number
  bdose <- readRDS(paste0("/auto/pmd-02/figi/HRC_BDose/FIGI_snpid_fix_chr", chr, ".rds"))
  dosages <- GetSNPValues(bdose, vector_of_snps, geneProb = F)
  out <- data.frame(dosages)
  #out$vcfid <- rownames(out)
  return(out)
}

out <- lapply(snps_list, GetSNPValues_wrapper)

out <- unname(out)
out <- do.call(cbind, out)
out$vcfid <- rownames(out)

saveRDS(out, file = paste0("FIGI_genotype_dosages_", filename, ".rds"))


