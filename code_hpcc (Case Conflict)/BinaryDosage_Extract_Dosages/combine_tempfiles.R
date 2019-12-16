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

# IOU - combine temporary files into a single data.frame for convenience. 

