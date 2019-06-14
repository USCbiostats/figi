#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti

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

#batch_list="corect_oncoarray"

for chr in {1..11}
do
    for b in ${batch_list}
    do
        plink --vcf ${b}_chr${chr}.recode.vcf --double-id --keep-allele-order --make-bed --out ./plink/${b}_chr${chr}
        cp ./plink/${b}_chr${chr}.bim ./plink/${b}_chr${chr}_backup.bim
        Rscript bimFix.R ./plink/${b}_chr${chr}
    done
done