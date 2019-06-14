#====================================================
# Testing the updated merge function
#
#====================================================

library(BinaryDosage)

#bd_file1 <- "/staging/dvc/andreeki/bdose/ukbiobank_chr5.bdose"
#bd_file1 <- "/auto/pmd-02/figi/HRC_BDose/ccfr_omni_chr5.bdose"
#bd_file2 <- "/auto/pmd-02/figi/HRC_BDose/mecc_chr5.bdose"

args <- commandArgs(trailingOnly=T)
chr <- args[1]


#etal <- list.files(path = "/auto/pmd-02/figi/HRC_BDose", pattern = "chr5", full.names = T)

etal <- c(paste0("/auto/pmd-02/figi/HRC_BDose/axiom_acs_aus_nf_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/axiom_mecc_cfr_ky_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/ccfr_1m_1mduo_reimpute_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/ccfr_omni_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/corect_oncoarray_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/mecc_chr", chr, ".bdose"))

etal <- c(paste0("/auto/pmd-02/figi/HRC_BDose/axiom_acs_aus_nf_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/axiom_mecc_cfr_ky_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/ccfr_1m_1mduo_reimpute_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/ccfr_omni_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/corect_oncoarray_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/corect_oncoarray_nonEUR_reimpute_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/mecc_chr", chr, ".bdose"))


#ukb <- "/staging/dvc/andreeki/bdose/ukbiobank_chr5.bdose"

#ukb_etal <- c(etal, ukb)

mergedFile <- paste0("/staging/dvc/andreeki/bdose/FIGI_CORECT_chr", chr, ".bdose")

MergeBD(mergedFile, etal)