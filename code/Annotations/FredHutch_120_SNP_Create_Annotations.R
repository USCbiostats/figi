#=======================================================================================#
# Process annotations 
# - previously reported GWAS
# - new hits reported by fred hutch
#
# Create files to use with EasyStrata
# Updated 06/03/2019
#=======================================================================================#

# single file, no LD based annotation, just top hits
annot <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_061918.tsv")[, -c('LL_95PCT_CI_BETA', 'UL_95PCT_CI_BETA')] %>% 
  filter(SUGGESTIVE == "no") %>% 
  dplyr::rename(Chr = CHROM,
                Pos = POSITION) %>% 
  mutate(Colour = "green", 
         SNP = paste(Chr, Pos, sep = ":"),
         ID = gsub("/", "_", VARIANT))


# i think stephanie sent old one by accident, need to ask her tomorrow. 



