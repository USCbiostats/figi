#====================================================
# Testing the updated merge function
#
#====================================================

library(BinaryDosage)

args <- commandArgs(trailingOnly=T)
chr <- args[1]

etal <- c(paste0("/auto/pmd-02/figi/HRC_BDose/axiom_acs_aus_nf_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/axiom_mecc_cfr_ky_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/ccfr_1m_1mduo_reimpute_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/ccfr_omni_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/corect_oncoarray_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/corect_oncoarray_nonEUR_reimpute_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/mecc_chr", chr, ".bdose"))


etal <- c(paste0("/auto/pmd-02/figi/HRC_BDose/corsa_axiom_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/cytosnp_comb_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/dachs3_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/initial_comb_datasets_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/newfoundland_omniquad_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/omni_comb_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/omniexpress_exomechip_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/oncoarray_to_usc_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/plco_3_chr", chr, ".bdose"),
          paste0("/auto/pmd-02/figi/HRC_BDose/reach_chr", chr, ".bdose"))



#ukb <- "/staging/dvc/andreeki/bdose/ukbiobank_chr5.bdose"

#ukb_etal <- c(etal, ukb)

mergedFile <- paste0("/auto/pmd-02/figi/HRC_BDose/FIGI_GECCO_chr", chr, ".bdose")

MergeBD(mergedFile, etal)
