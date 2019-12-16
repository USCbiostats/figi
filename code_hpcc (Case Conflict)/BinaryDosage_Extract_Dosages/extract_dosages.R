#=============================================================================#
# FIGI - Extract Dosages from BinaryDosage files
# 11/03/2019
#
# Take a long vector of index positions (CHROMOSOME SPECIFIC)
#   split vector into chunks, then loop (or apply) over them to extract dosages
#   output chunked dosages into rds files
#   combine rds files into a single data.frame for convenience
#=============================================================================#
library(BinaryDosage)

# always take two arguments - chr and snplist (vector of index positions)
args <- commandArgs(trailingOnly=T)
chr <- args[1] # chromosome 
filename <- args[2] # just a string of the file name containing SNP index positions in BinaryDosage files

snpsToRead <- readRDS(paste0(filename, ".rds"))
snpsToRead_chunk <- split(snpsToRead, ceiling(seq_along(snpsToRead)/500))


# if things start crashing, might have to resort to creating multiple slurm jobs..
bdose <- readRDS(paste0("/auto/pmd-02/figi/HRC_BDose/FIGI_chr", chr, ".rds")) # FIGI_chr.rds index file

#GetSNPValues_wrapper <- function(chunk) {
#  tmp <- GetSNPValues(bdose, chunk, geneProb = F)
#  out <- data.frame(tmp)
#  out$vcfid <- rownames(out)
#  saveRDS(out, file = paste0(filename, "_TEMP_", chunk, ".rds"))
#}

#for(x in snpsToRead_chunk) {
#  GetSNPValues_wrapper(x)
#}

counter <- 1
for(x in snpsToRead_chunk) {
    tmp <- GetSNPValues(bdose, x, geneProb = F)
    out <- data.frame(tmp)
    out$vcfid <- rownames(out)
    saveRDS(out, file = paste0(filename, "_TMP", counter, ".rds"))
    counter <- counter + 1
}

