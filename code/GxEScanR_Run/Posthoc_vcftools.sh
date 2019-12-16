#!/bin/bash

# very low P values, investigate
#22:17265194 

OUT=/staging/dvc/andreeki/posthoc
REF=/auto/pmd-02/figi/HRC_VCF_SampleRename

cd ${OUT}

wtf_gxescan () {
echo "#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti

vcftools --gzvcf ${REF}/$1_chr$2.vcf.gz --positions ~/FIGI_ALL_chr22_topP.txt --recode --out ${OUT}/$1_topSNPs

" | sbatch
}

for batch in \
"axiom_acs_aus_nf" \
"axiom_mecc_cfr_ky" \
"ccfr_1m_1mduo" \
"ccfr_omni" \
"corect_oncoarray" \
"corect_oncoarray_nonEUR" \
"corsa_axiom" \
"cytosnp_comb" \
"dachs3" \
"initial_comb_datasets" \
"mecc" \
"newfoundland_omniquad" \
"omni_comb" \
"omniexpress_exomechip" \
"oncoarray_to_usc" \
"plco_3" \
"reach" \
"ukbiobank" \   
do
    wtf_gxescan $batch 22
done




for file in /dir/*
do
  bgzip -c ${file} > ${file}.gz
  tabix ${file}.gz
done



# next time just use bcftools for this kind of thing. vcftools is so clunky

vcf-merge \
axiom_acs_aus_nf_topSNPs.recode.vcf.gz \
corsa_axiom_topSNPs.recode.vcf.gz \
omni_comb_topSNPs.recode.vcf.gz \
axiom_mecc_cfr_ky_topSNPs.recode.vcf.gz \
cytosnp_comb_topSNPs.recode.vcf.gz \
omniexpress_exomechip_topSNPs.recode.vcf.gz \
ccfr_1m_1mduo_topSNPs.recode.vcf.gz \
dachs3_topSNPs.recode.vcf.gz \
oncoarray_to_usc_topSNPs.recode.vcf.gz \
ccfr_omni_topSNPs.recode.vcf.gz \
initial_comb_datasets_topSNPs.recode.vcf.gz  \
plco_3_topSNPs.recode.vcf.gz \
corect_oncoarray_nonEUR_topSNPs.recode.vcf.gz  \
mecc_topSNPs.recode.vcf.gz  \
reach_topSNPs.recode.vcf.gz \
corect_oncoarray_topSNPs.recode.vcf.gz \
newfoundland_omniquad_topSNPs.recode.vcf.gz  \
ukbiobank_topSNPs.recode.vcf.gz > FIGI_ALL_chr22_topSNPs.vcf

bgzip -c FIGI_ALL_chr22_topSNPs.vcf > FIGI_ALL_chr22_topSNPs.vcf.gz
tabix FIGI_ALL_chr22_topSNPs.vcf.gz