#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --output=./logs/j03_mergeALL.log

REF=/auto/pmd-02/figi/PCA
OUT=/staging/dvc/andreeki/pca_ibd
cd ${OUT}

# Change SNP names for ukbiobank and reach)
# (Already performed in previous PCA calculation runs)

# cp ukbiobank_backbone.bim ukbiobank_backbone.bim.bck
# awk '{print $1, $1":"$4, $3, $4, $5, $6}' OFS='\t' ukbiobank_backbone.bim.bck > ukbiobank_backbone.bim
# rm ukbiobank_backbone.bim.bck

# cp reach_backbone.bim reach_backbone.bim.bck
# awk '{print $1, $1":"$4, $3, $4, $5, $6}' OFS='\t' reach_backbone.bim.bck > reach_backbone.bim
# rm reach_backbone.bim.bck


# Create mergelist to use in plink command
batch_list="axiom_acs_aus_nf_backbone \
axiom_mecc_cfr_ky_backbone \
ccfr_1m_1mduo_reimpute_backbone \
ccfr_omni_backbone \
corect_oncoarray_backbone \
corect_oncoarray_nonEUR_reimpute_backbone \
corsa_axiom_backbone \
cytosnp_comb_backbone \
dachs3_backbone \
initial_comb_datasets_backbone \
mecc_backbone \
newfoundland_omniquad_backbone \
omni_comb_backbone \
omniexpress_exomechip_backbone \
oncoarray_to_usc_backbone \
plco_3_backbone \
reach_backbone \
ukbiobank_backbone"

for batch in ${batch_list}
do 
    echo ${REF}/${batch}
done > ${OUT}/tmp/mergelist_all.txt

plink --merge-list ${OUT}/tmp/mergelist_all.txt --make-bed --out ${OUT}/TMP1

#------ Important Note ------#
# there's less than 100% genotyping rate, that's because of oncoarray_to_usc incomplete samples across all chromosome files. Just remember that when you subset theh GWAS/GxE sets, then genotyping rates becomes 100%. 

#------ IMPORTANT NOTES ------#
# these markers on ukbiobank imputation batch:
#12:99935010
#6:41924033

# remove these from chromosome specific plink file (during chromosome merging e.g. before running j02.sh):
#6:41924033_A_AAAAAT
#12:99935010_AAAAC_A

# (i'm going to remove using command line because i had to use plink2 to convert ukbiobank, and this version doesn't allow --snp-only option yet)
# --- meaning - removal of these variants isn't documented here but trust the fact that they're removed from the merged plink files


# important - oncoarray_to_usc... not all chromosomes have exactly same samples. 
# usually kept only intersecting samples
# some might not be found in samplefile
# Decision - keep only samples listed in the samplefile
# again, used command line ("oncoarray_to_usc_samplekeep.txt")

# Create GWAS and GxE Sets
plink --bfile ${OUT}/TMP1 --keep FIGI_PCA_GwasSet.txt --make-bed --out ${OUT}/FIGI_GwasSet
plink --bfile ${OUT}/TMP1 --keep FIGI_PCA_GxESet.txt --make-bed --out ${OUT}/FIGI_GxESet


#plink --bfile TMP1 --keep FIGI_PCA_GwasSet.txt --make-bed --out FIGI_GwasSet
#plink --bfile TMP1 --keep FIGI_PCA_GxESet.txt --make-bed --out FIGI_GxESet
