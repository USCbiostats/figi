#!/bin/bash

batch_list="axiom_acs_aus_nf \
axiom_mecc_cfr_ky \
ccfr_1m_1mduo_reimpute \
ccfr_omni \
corect_oncoarray \
corect_oncoarray_nonEUR_reimpute \
corsa_axiom \
cytosnp_comb \
dachs3 \
initial_comb_datasets \
mecc \
newfoundland_omniquad \
omni_comb \
omniexpress_exomechip \
oncoarray_to_usc \
plco_3 \
reach \
ukbiobank"

batch_list="axiom_acs_aus_nf \
axiom_mecc_cfr_ky \
ccfr_1m_1mduo_reimpute \
corect_oncoarray \
cytosnp_comb \
dachs3 \
initial_comb_datasets \
mecc \
omni_comb \
omniexpress_exomechip \
oncoarray_to_usc \
plco_3 \
reach"

#batch_list="axiom_mecc_cfr_ky"

for chr in {20..20}
do
    for b in ${batch_list}
    do
        #chr=21
        sbatch --job-name=${b}_chr${chr}_plink_bed.run --output=./logs/${b}_chr${chr}_plink_bed.log --export=batch=${b},chr=${chr} 02_bed.sh
    done
done
#sbatch --export=batch='axiom_acs_aus_nf',chr=22 vcftools_pp.sh