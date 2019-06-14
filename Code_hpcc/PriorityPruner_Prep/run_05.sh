#!/bin/bash

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

#batch_list="axiom_acs_aus_nf"

### Create list of files to merge using plink (one time run only)
#for chr in {1..22}
#do
#    for batch in ${batch_list}
#    do
#        echo ./plink/${batch}_chr${chr}_filter
#    done > ./plink/merge_list_chr${chr}.txt
#done

for chr in {1..22}
do
    sbatch --job-name=pp_chr${chr}.run --output=./logs/pp_chr${chr}.log --export=chr=${chr} 05_priority_pruner.sh
done