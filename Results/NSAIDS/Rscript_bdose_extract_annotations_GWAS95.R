annot <- fread("~/bin/EasyStrata/crc_gwas_125k_indep_signals_061918.tsv")[, -c('LL_95PCT_CI_BETA', 'UL_95PCT_CI_BETA')] %>% 
  filter(SUGGESTIVE == "no") %>% 
  dplyr::rename(Chr = CHROM,
                Pos = POSITION) %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"))


for(chr in 1:22) {
  snpsToGet <- annot[which(annot$Chr == chr), 'SNP']
  saveRDS(snpsToGet, file = paste0("~/Dropbox/code/FIGI_GxEScanR_Simple/working/GWAS_95loci_chr", chr, ".rds"))
}

snpsToGet <- annot[which(annot$Chr == 20), 'SNP']
