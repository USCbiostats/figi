#=============================================================================#
# After calculating Rsq using the forked BinaryDosage package (see my github)
#
# Need to process matrices:
# 1) match index position using BDose index files to get SNP ID (with alt/ref info)
# 2) create data.frames for each chromosome
# 3) filter now, or just filter before generating results plots etc
#=============================================================================#

# start with chromosome 22
index <- readRDS("~/data/FIGI_BDose_IndexFiles/FIGI_chr22.rds")[[11]] %>% 
  dplyr::mutate(SNP = paste(SNPID, Reference, Alternate, sep = ":"),
                Rsq_AA = readRDS("~/data/Rsq_Estimate/FIGI_RsqEstimate_chr22.rds"))


# genome-wide
index_function <- function(chr) {
  readRDS(paste0("~/data/FIGI_BDose_IndexFiles/FIGI_chr", chr, ".rds"))[[11]] %>% 
    dplyr::mutate(SNP = paste(SNPID, Reference, Alternate, sep = ":"),
                  Rsq_AA = readRDS(paste0("~/data/Rsq_Estimate/FIGI_RsqEstimate_chr", chr, ".rds")))
}

rsq_aa_genomewide <- do.call(rbind, lapply(as.list(1:22), index_function))
saveRDS(rsq_aa_genomewide, file = "/media/work/data/Rsq_Estimates_AA/rsq_aa_genomewide.rds", version = 2)

rsq_aa_genomwide_filter <- rsq_aa_genomewide %>% 
  dplyr::filter(Rsq_AA >= 0.8)
saveRDS(rsq_aa_genomwide_filter, file = "/media/work/data/Rsq_Estimates_AA/rsq_aa_genomewide_filter_pt8.rds", version = 2)
