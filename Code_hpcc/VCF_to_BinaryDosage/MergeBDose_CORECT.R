#====================================================
# Merge BinaryDosage Files
# 
# CORECT study imputation batches:
# axiom_acs_aus_nf
# axiom_mecc_cfr_ky
# ccfr_1m_1mduo_reimpute
# ccfr_omni
# corect_oncoarray
# {corect_oncoarray_nonEUR_reimpute .. this imputation has issues maybe exclude}
# mecc
#
#====================================================
library(BinaryDosage)

# command argument - chromosome
args <- commandArgs(trailingOnly=T)
chr <- args[1]

batchlist <- c(paste0("/auto/pmd-02/figi/HRC_BDose/axiom_acs_aus_nf_chr", chr, ".bdose"),
               paste0("/auto/pmd-02/figi/HRC_BDose/axiom_mecc_cfr_ky_chr", chr, ".bdose"),
               paste0("/auto/pmd-02/figi/HRC_BDose/ccfr_1m_1mduo_reimpute_chr", chr, ".bdose"),
               paste0("/auto/pmd-02/figi/HRC_BDose/ccfr_omni_chr", chr, ".bdose"),
               paste0("/auto/pmd-02/figi/HRC_BDose/corect_oncoarray_chr", chr, ".bdose"),
               paste0("/auto/pmd-02/figi/HRC_BDose/mecc_chr", chr, ".bdose"))

#etal <- c(paste0("/auto/pmd-02/figi/HRC_BDose/axiom_acs_aus_nf_chr", chr, ".bdose"),
#          paste0("/auto/pmd-02/figi/HRC_BDose/axiom_mecc_cfr_ky_chr", chr, ".bdose"),
#          paste0("/auto/pmd-02/figi/HRC_BDose/ccfr_1m_1mduo_reimpute_chr", chr, ".bdose"),
#          paste0("/auto/pmd-02/figi/HRC_BDose/ccfr_omni_chr", chr, ".bdose"),
#          paste0("/auto/pmd-02/figi/HRC_BDose/corect_oncoarray_chr", chr, ".bdose"),
#          paste0("/auto/pmd-02/figi/HRC_BDose/corect_oncoarray_nonEUR_reimpute_chr", chr, ".bdose"),
#          paste0("/auto/pmd-02/figi/HRC_BDose/mecc_chr", chr, ".bdose"))

mergedFile <- paste0("/auto/pmd-02/figi/HRC_BDose/FIGI_CORECT_chr", chr, ".bdose")
MergeBD(mergedFile, batchlist)
