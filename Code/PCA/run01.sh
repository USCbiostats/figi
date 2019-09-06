#!/bin/bash


#sbatch --output=./logs/j01_axiom_acs_aus_nf_chr%a.log --export=batch=axiom_acs_aus_nf j01.sh
#sbatch --output=./logs/j01_axiom_mecc_cfr_ky_chr%a.log --export=batch=axiom_mecc_cfr_ky j01.sh
#sbatch --output=./logs/j01_ccfr_1m_1mduo_reimpute_chr%a.log --export=batch=ccfr_1m_1mduo_reimpute j01.sh
#sbatch --output=./logs/j01_ccfr_omni_chr%a.log --export=batch=ccfr_omni j01.sh

sbatch --output=./logs/j01_corect_oncoarray_chr%a.log --export=batch=corect_oncoarray j01.sh
sbatch --output=./logs/j01_corect_oncoarray_nonEUR_reimpute_chr%a.log --export=batch=corect_oncoarray_nonEUR_reimpute j01.sh
sbatch --output=./logs/j01_corsa_axiom_chr%a.log --export=batch=corsa_axiom j01.sh
sbatch --output=./logs/j01_cytosnp_comb_chr%a.log --export=batch=cytosnp_comb j01.sh


#sbatch --output=./logs/j01_dachs3_chr%a.log --export=batch=dachs3 j01.sh
#sbatch --output=./logs/j01_initial_comb_datasets_chr%a.log --export=batch=initial_comb_datasets j01.sh
#sbatch --output=./logs/j01_mecc_chr%a.log --export=batch=mecc j01.sh
#sbatch --output=./logs/j01_newfoundland_omniquad_chr%a.log --export=batch=newfoundland_omniquad j01.sh


#sbatch --output=./logs/j01_omni_comb_chr%a.log --export=batch=omni_comb j01.sh
#sbatch --output=./logs/j01_omniexpress_exomechip_chr%a.log --export=batch=omniexpress_exomechip j01.sh
#sbatch --output=./logs/j01_oncoarray_to_usc_chr%a.log --export=batch=oncoarray_to_usc j01.sh
#sbatch --output=./logs/j01_plco_3_chr%a.log --export=batch=plco_3 j01.sh

#sbatch --output=./logs/j01_reach_chr%a.log --export=batch=reach j01.sh
#sbatch --output=./logs/j01_UKB1_chr%a.log --export=batch=UKB1 j01.sh
#sbatch --output=./logs/j01_UKB2_chr%a.log --export=batch=UKB2 j01.sh
