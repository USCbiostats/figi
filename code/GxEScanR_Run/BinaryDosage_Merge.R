library(BinaryDosage)
library(GxEScanR)

args <- commandArgs(trailingOnly=T)
chr <- args[1]

setwd("/auto/pmd-02/figi/HRC_BDose")

# files <-  c(paste0("axiom_acs_aus_nf_chr", chr, ".bdose"),
#             paste0("axiom_mecc_cfr_ky_chr", chr, ".bdose"),
#             paste0("ccfr_1m_1mduo_reimpute_chr", chr, ".bdose"),
#             paste0("ccfr_omni_chr", chr, ".bdose"),
#             paste0("corect_oncoarray_chr", chr, ".bdose"),
#             paste0("corect_oncoarray_nonEUR_reimpute_chr", chr, ".bdose"),
#             paste0("corsa_axiom_chr", chr, ".bdose"),
#             paste0("cytosnp_comb_chr", chr, ".bdose"),
#             paste0("dachs3_chr", chr, ".bdose"),
#             paste0("initial_comb_datasets_chr", chr, ".bdose"),
#             paste0("mecc_chr", chr, ".bdose"),
#             paste0("newfoundland_omniquad_chr", chr, ".bdose"),
#             paste0("omni_comb_chr", chr, ".bdose"),
#             paste0("omniexpress_exomechip_chr", chr, ".bdose"),
#             paste0("oncoarray_to_usc_chr", chr, ".bdose"),
#             paste0("plco_3_chr", chr, ".bdose"),
#             paste0("reach_chr", chr, ".bdose"),
#             paste0("ukbiobank_chr", chr, ".bdose"))

files <-  c(paste0("axiom_acs_aus_nf_chr", chr, ".bdose"),
            paste0("axiom_mecc_cfr_ky_chr", chr, ".bdose"),
            paste0("ccfr_1m_1mduo_reimpute_chr", chr, ".bdose"),
            paste0("ccfr_omni_chr", chr, ".bdose"),
            paste0("corect_oncoarray_chr", chr, ".bdose"),
            paste0("corect_oncoarray_nonEUR_reimpute_chr", chr, ".bdose"),
            paste0("corsa_axiom_chr", chr, ".bdose"),
            paste0("cytosnp_comb_chr", chr, ".bdose"),
            paste0("dachs3_chr", chr, ".bdose"),
            paste0("initial_comb_datasets_chr", chr, ".bdose"),
            paste0("mecc_chr", chr, ".bdose"),
            paste0("newfoundland_omniquad_chr", chr, ".bdose"),
            paste0("omni_comb_chr", chr, ".bdose"),
            paste0("omniexpress_exomechip_chr", chr, ".bdose"),
            paste0("oncoarray_to_usc_chr", chr, ".bdose"),
            paste0("plco_3_chr", chr, ".bdose"),
            paste0("reach_chr", chr, ".bdose"))

start_time <- Sys.time()
MergeBD42(paste0("/staging/dvc/andreeki/figi/FIGI_ALL_chr", chr, ".bdose"), files)
end_time <- Sys.time()

start_time
end_time
end_time - start_time
