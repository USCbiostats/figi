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
    plink --merge-list ./plink/merge_list_chr${chr}.txt --keep-allele-order --make-bed --out ./plink/final_chr${chr}
done


#vcftools --gzvcf /auto/pmd-02/figi/HRC_VCF_SampleRename/${batch}_chr${chr}.vcf.gz --positions PriorityPruner_noUKB_chr${chr}_snplist.txt --keep PriorityPruner_noUKB_10000_samplelist.txt --recode --out ${batch}_chr${chr}