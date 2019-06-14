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
        plink --bfile ./plink/${b}_chr${chr} --extract ./plink/PriorityPruner_noUKB_ALL_snplist.txt  --make-bed --out ./plink/${b}_chr${chr}_filter
    done
done


#vcftools --gzvcf /auto/pmd-02/figi/HRC_VCF_SampleRename/${batch}_chr${chr}.vcf.gz --positions PriorityPruner_noUKB_chr${chr}_snplist.txt --keep PriorityPruner_noUKB_10000_samplelist.txt --recode --out ${batch}_chr${chr}