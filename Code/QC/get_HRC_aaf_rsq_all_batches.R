#=============================================================================#
# Obtain and organize Rsq and MAF info from HRC Info files
# in order to incorporate quality filter for GxEScan results
# 
# some notes:
# - dachs is part of another imputation batch
# - minimac4 reports different snp IDs than minimac3, process separately
#   (this will not be an issue in a few weeks)
#   (but be mind that you might need to redo these steps to generate avg rsq again)
#=============================================================================#
library(tidyverse)
library(data.table)

batch_list_minimac3 <-  list('axiom_acs_aus_nf',
                             'axiom_mecc_cfr_ky',
                             'ccfr_1m_1mduo_reimpute',
                             'ccfr_omni',
                             'corect_oncoarray',
                             'corsa_axiom',
                             'cytosnp_comb',
                             'initial_comb_datasets',
                             'mecc',
                             'newfoundland_omniquad',
                             'omni_comb',
                             'omniexpress_exomechip',
                             'oncoarray_to_usc',
                             'plco_3')

batch_list_minimac4 <- list('reach', 'ukbiobank')

# only difference is recreating the ID field
wrap_minimac3 <- function(x, chr){
  varname1 <- paste0(x, "_ALT_Frq")
  varname2 <- paste0(x, "_Rsq")
  varname3 <- paste0(x, "_ID")
  
  # for future ref, this is called dynamic var names
  fread(paste0("/home/rak/data/HRC_InfoFile/", x, "/chr", chr, ".info.gz")) %>% 
    mutate(ID = paste(SNP, `REF(0)`, `ALT(1)`, sep = ":"),
           !!varname1 := ALT_Frq,
           !!varname2 := Rsq) %>% 
    # filter(!!varname3 %in% snplist) %>% 
    dplyr::select(ID, !!varname1, !!varname2)
}


wrap_minimac4 <- function(x, chr){
  varname1 <- paste0(x, "_ALT_Frq")
  varname2 <- paste0(x, "_Rsq")
  varname3 <- paste0(x, "_ID")
  
  fread(paste0("/home/rak/data/HRC_InfoFile/", x, "/chr", chr, ".info.gz")) %>% 
    mutate(ID = SNP,
           !!varname1 := ALT_Frq,
           !!varname2 := Rsq) %>% 
    # filter(!!varname3 %in% snplist) %>% 
    dplyr::select(ID, !!varname1, !!varname2)
}


# reduce does it sequentially (that's what you want)
# Save files in ~/data/HRC_InfoFile_Merged
for(chr in c(22:22)) {
  test <- Reduce(inner_join, lapply(batch_list1, wrap_minimac3, chr = chr))
  test2 <- Reduce(inner_join, lapply(batch_list2, wrap_minimac4, chr = chr))
  z <- inner_join(test, test2, by = "ID")
  saveRDS(z, file = paste0("~/data/HRC_InfoFile_Merged/HRC_InfoFile_Chr", chr, ".rds"))
}


#-----------------------------------------------------------------------------#
# calculate weighted avg Rsq across 16 batches
# save object as a vector of SNP IDs of format CHR:POS:REF:ALT
#-----------------------------------------------------------------------------#

batch_list <- c("axiom_acs_aus_nf",
                "axiom_mecc_cfr_ky",
                "ccfr_1m_1mduo_reimpute",
                "ccfr_omni",
                "corect_oncoarray",
                "corsa_axiom",
                "cytosnp_comb",
                "initial_comb_datasets",
                "mecc",
                "newfoundland_omniquad",
                "omni_comb",
                "omniexpress_exomechip",
                "oncoarray_to_usc",
                "plco_3",
                "reach",
                "ukbiobank")
sample_sizes <- c(2766, 7501, 2180, 1600, 36621, 2467, 10908, 5924, 983, 637, 5986, 5503, 20912, 4864, 750, 27594)
chroms <- names(table(results_GxE$Chromosome))

# get info from each chromosome info file
# calculate by chromosome
chr <- 22
avg_rsq_chrom <- function(chr) {
    tmp1 <- readRDS(paste0("~/data/HRC_InfoFile_Merged/HRC_InfoFile_Chr", chr, ".rds"))
    idlist <- tmp1[, c("ID")]
    tmp2 <- tmp1 %>% dplyr::select(contains("Rsq")) %>% 
      mutate_if(is.character, as.numeric)
    tmp3 <- data.frame(mapply(`*`, tmp2, sample_sizes)) %>% 
      mutate(Rsq_tot = rowSums(.),
             Rsq_avg = Rsq_tot / sum(sample_sizes),
             ID = idlist) %>% 
      dplyr::select(ID, Rsq_avg)
}

# apply function over list of chromosomes, then just concatenate
hrc_Rsq_avg <- do.call(rbind, lapply(as.list(1:22), avg_rsq_chrom))

# maybe filter out some very rare alleles from the above... 
hrc_MAF <- readRDS("~/data/HRC/HRC.r1-1.rds")
hrc_MAF_filter <- hrc_MAF %>% 
  filter(maf > 0.005) %>% 
  mutate(ID = paste(V1, V2, V3, V4, sep = ":"))

x <- inner_join(hrc_MAF_filter[, c("ID", "maf", "mafbin")], hrc_Rsq_avg, by = 'ID')

saveRDS(x, file = "~/data/HRC_InfoFile_Merged/HRC_WtAvgRsq_HRCRefPanelMAFs.rds", version = 2)


wt <- readRDS("~/data/HRC_InfoFile_Merged/HRC_InfoFile_Chr22.rds")
