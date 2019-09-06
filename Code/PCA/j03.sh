#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB
#SBATCH --output=./logs/j03_mergeALL.log

DIR=/auto/pmd-02/figi/PCA
cd ${OUT}

# Change SNP names for ukbiobank and reach)
# (old - no longer necessary)

# cp ukbiobank_backbone.bim ukbiobank_backbone.bim.bck
# awk '{print $1, $1":"$4, $3, $4, $5, $6}' OFS='\t' ukbiobank_backbone.bim.bck > ukbiobank_backbone.bim
# rm ukbiobank_backbone.bim.bck

# cp reach_backbone.bim reach_backbone.bim.bck
# awk '{print $1, $1":"$4, $3, $4, $5, $6}' OFS='\t' reach_backbone.bim.bck > reach_backbone.bim
# rm reach_backbone.bim.bck


# Create mergelist to use in plink command
#batch_list="axiom_acs_aus_nf_backbone \
#axiom_mecc_cfr_ky_backbone \
#ccfr_1m_1mduo_reimpute_backbone \
#ccfr_omni_backbone \
#corect_oncoarray_backbone \
#corect_oncoarray_nonEUR_reimpute_backbone \
#corsa_axiom_backbone \
#cytosnp_comb_backbone \
#dachs3_backbone \
#initial_comb_datasets_backbone \
#mecc_backbone \
#newfoundland_omniquad_backbone \
#omni_comb_backbone \
#omniexpress_exomechip_backbone \
#oncoarray_to_usc_backbone \
#plco_3_backbone \
#reach_backbone \
#ukbiobank_backbone"

#for batch in ${batch_list}
#do 
#    echo ${REF}/${batch}
#done > ${OUT}/files/mergelist_all.txt

echo  "/auto/pmd-02/figi/PCA/files/axiom_acs_aus_nf_backbone
/auto/pmd-02/figi/PCA/files/axiom_mecc_cfr_ky_backbone
/auto/pmd-02/figi/PCA/files/ccfr_1m_1mduo_reimpute_backbone
/auto/pmd-02/figi/PCA/files/ccfr_omni_backbone
/auto/pmd-02/figi/PCA/files/corect_oncoarray_backbone
/auto/pmd-02/figi/PCA/files/corect_oncoarray_nonEUR_reimpute_backbone
/auto/pmd-02/figi/PCA/files/corsa_axiom_backbone
/auto/pmd-02/figi/PCA/files/cytosnp_comb_backbone
/auto/pmd-02/figi/PCA/files/dachs3_backbone
/auto/pmd-02/figi/PCA/files/initial_comb_datasets_backbone
/auto/pmd-02/figi/PCA/files/mecc_backbone
/auto/pmd-02/figi/PCA/files/newfoundland_omniquad_backbone
/auto/pmd-02/figi/PCA/files/omni_comb_backbone
/auto/pmd-02/figi/PCA/files/omniexpress_exomechip_backbone
/auto/pmd-02/figi/PCA/files/oncoarray_to_usc_backbone
/auto/pmd-02/figi/PCA/files/plco_3_backbone
/auto/pmd-02/figi/PCA/files/reach_backbone
/auto/pmd-02/figi/PCA/files/UKB1_backbone
/auto/pmd-02/figi/PCA/files/UKB2_backbone" > ${DIR}/j03_mergelist_all.txt

plink --merge-list ${DIR}/j03_mergelist_all.txt --make-bed --out ${DIR}/TMP1

#------ Important Note ------#
#(no longer applies, removed sample inconsistency in oncoarray_to_usc))
# there's less than 100% genotyping rate, that's because of oncoarray_to_usc incomplete samples across all chromosome files. Just remember that when you subset theh GWAS/GxE sets, then genotyping rates becomes 100%. 

#------ IMPORTANT NOTES ------#
#(no longer applies - UKB was reimputed with minimac3)
# these markers on ukbiobank imputation batch:
#12:99935010
#6:41924033

# remove these from chromosome specific plink file (during chromosome merging e.g. before running j02.sh):
#6:41924033_A_AAAAAT
#12:99935010_AAAAC_A

# (i'm going to remove using command line because i had to use plink2 to convert ukbiobank, and this version doesn't allow --snp-only option yet)
# --- meaning - removal of these variants isn't documented here but trust the fact that they're removed from the merged plink files

#(also not applicable, see above)
# important - oncoarray_to_usc... not all chromosomes have exactly same samples. 
# usually kept only intersecting samples
# some might not be found in samplefile
# Decision - keep only samples listed in the samplefile
# again, used command line ("oncoarray_to_usc_samplekeep.txt")

# Create GWAS and GxE Sets
plink --bfile ${DIR}/TMP1 --keep FIGI_PCA_GwasSet_190729.txt --make-bed --out ${DIR}/FIGI_GwasSet_190729
plink --bfile ${DIR}/TMP1 --keep FIGI_PCA_GxESet_190729.txt --make-bed --out ${DIR}/FIGI_GxESet_190729


