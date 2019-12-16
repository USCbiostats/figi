#=======================================================================================#
# Process annotations 
# - previously reported GWAS
# - new hits reported by fred hutch
#
# Note: keep in mind that allele filters depend on exposure and covariates
# so some previously reported GWAS hits and/or newer findings from FredHutch 
# might be missing from your own analysis... 
#
#
# Create files to use with EasyStrata
#=======================================================================================#

# two sets of annotations - previously reported in literature vs newly reported in Huyghe 2019
jh_lit <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_061918.tsv")[, -c('LL_95PCT_CI_BETA', 'UL_95PCT_CI_BETA')] %>% 
  filter(SUGGESTIVE == "no",
         PREVIOUSLY_REPORTED == "yes") %>% 
  dplyr::rename(Chr = CHROM,
                Pos = POSITION) %>% 
  mutate(Colour = "green", 
         SNP = paste(Chr, Pos, sep = ":"),
         ID = gsub("/", "_", VARIANT))

# write.table(jh_lit, file = "~/data/Annotations/crc_gwas_125k_indep_signals_95_GWAS.tsv", quote = F, row.names = F, sep = '\t')
# check_lit <- anti_join(jh_lit, results_G, by = 'ID') # likely filtered out of our results during GxEScanR.. missing some important ones?
# filter(results_G, Location %in% check$Location)
# filter(ukb_filter, Location %in% check$Pos)


jh_new <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_061918.tsv")[, -c('LL_95PCT_CI_BETA', 'UL_95PCT_CI_BETA')] %>% 
  filter(SUGGESTIVE == "no",
         PREVIOUSLY_REPORTED == "no") %>% 
  dplyr::rename(Chr = CHROM,
                Pos = POSITION) %>% 
  mutate(Colour = "red", 
         SNP = paste(Chr, Pos, sep = ":"),
         ID = gsub("/", "_", VARIANT))

# write.table(jh_new, file = "~/data/Annotations/crc_gwas_125k_indep_signals_95_NEW.tsv", quote = F, row.names = F, sep = '\t')
# check <- anti_join(jh_new, results_G, by = 'ID') # filtered out too? How common/rare are these.. 
# filter(results_G, Location %in% check$Location) 
# filter(ukb_filter, Location %in% check$Pos)


jh <- rbind(jh_lit, jh_new)
write.table(jh, file = "~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata.tsv", quote = F, row.names = F, sep = '\t')


#-----------------------------------------------------------------------------#
# LD Based Annotations
#-----------------------------------------------------------------------------#

gwas95 <- readRDS("~/data/Annotations/Known_95_Loci_NatGenList_Jan2019.rds")

# quick test
# test <- gwas95[[30]] %>% 
#   dplyr::rename(Pos = end) %>% 
#   dplyr::mutate(Chr = gsub("chr", "", chr),
#                 Colour = case_when(0 <= r2 & r2 < 0.2 ~ "navy",
#                             0.2 <= r2 & r2 < 0.4 ~ "lightblue", 
#                             0.4 <= r2 & r2 < 0.6 ~ "green", 
#                             0.6 <= r2 & r2 < 0.8 ~ "orange", 
#                             0.8 <= r2 & r2 <= 1.0 ~ "red")) %>% 
#   dplyr::select(Chr, Pos, Colour)
# table(test$Colour)
# write.table(test, file = "~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv", quote = F, row.names = F, sep = '\t')

wtf <- do.call(rbind, lapply(gwas95, nrow))
wtf2 <- rep(names(gwas95), wtf)

# now do this for every marker region
gwas95_df <- do.call(rbind, gwas95) %>% 
  dplyr::rename(Pos = end) %>% 
  dplyr::mutate(Chr = gsub("chr", "", chr),
                Colour = case_when(0 <= r2 & r2 < 0.2 ~ "lightblue",
                                   0.2 <= r2 & r2 < 0.4 ~ "lightblue", 
                                   0.4 <= r2 & r2 < 0.6 ~ "green", 
                                   0.6 <= r2 & r2 < 0.8 ~ "orange", 
                                   0.8 <= r2 & r2 <= 1.0 ~ "red"),
                marker = wtf2)
  # dplyr::select(Chr, Pos, Colour, marker)
table(gwas95_df$Colour)
write.table(gwas95_df, file = "~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv", quote = F, row.names = F, sep = '\t')



# let's explore overlaps with my results (asp_ref)
gwas95_original <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_061918.tsv")[, -c('LL_95PCT_CI_BETA', 'UL_95PCT_CI_BETA')] %>% 
  filter(SUGGESTIVE == "no")


gwas95_df <- do.call(rbind, gwas95) %>% 
  filter(r2 == 1)
