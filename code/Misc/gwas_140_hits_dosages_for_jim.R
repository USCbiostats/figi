# gwas main hits csv file for jim

list.files("/media/work/FIGI/gwas_140_hits_dosages/")
x <- lapply(list.files("/media/work/FIGI/gwas_140_hits_dosages/", full.names = T, pattern = "TMP1"), readRDS)
y <- Reduce(inner_join, x)

write.csv(y, file = "/media/work/FIGI/gwas_140_hits_dosages/figi_gwas_140_hits_dosages.csv", quote = F, row.names = F)
