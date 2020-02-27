source("./multiallelic_helper_gwas.R")

args <- commandArgs(trailingOnly=T)
chr <- args[1]
exposure <- args[2]


wrap(chr, exposure)
