#====================================================
# FIGI GxEScan run on Merged Files
#
# John test - GetVCFInfo
#====================================================
library(BinaryDosage)

# args
args <- commandArgs(trailingOnly=T)
BDoseFile <- args[1]

# run GxEScan
start_time <- Sys.time()
bdose <- GetBDoseInfo(paste0(BDoseFile, ".bdose"), index = T)
end_time <- Sys.time()

saveRDS(bdose, file = paste0(BDoseFile, ".rds"))

message(paste("GetBDoseInfo runtime: ", difftime(end_time, start_time, units = 'mins')))
