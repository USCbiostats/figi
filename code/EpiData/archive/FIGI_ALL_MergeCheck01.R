library(GxEScanR)

args <- commandArgs(trailingOnly=T)
ibatch <- args[1]
chr <- args[2]

ref <- "/auto/pmd-02/figi/HRC_BDose/"
out <- "/staging/dvc/andreeki/bdose/qc/"
filename <- paste0(ibatch, "_chr", chr)

bdose <- GxEScanR::GetBinaryDosageInfo(paste0(ref, filename, ".bdose"))
write.table(bdose[[7]], file = paste0(out, filename, "_samplelist.txt"), quote = F, row.names = F)