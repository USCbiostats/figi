#!/bin/bash
# merge chromosomes

OUT=/staging/dvc/andreeki/pca_ibd
cd ${OUT}

merge_backbone () {
echo "#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti

cd ${OUT}

plink --merge-list ${OUT}/mergelist_$1.txt --make-bed --out ${OUT}/$1_backbone

" | sbatch
}


# "axiom_acs_aus_nf"
# "axiom_mecc_cfr_ky"
# "ccfr_1m_1mduo"
# "ccfr_omni"
# "corect_oncoarray"
# "corect_oncoarray_nonEUR"
# "corsa_axiom"
# "cytosnp_comb"
# "dachs3"
# "initial_comb_datasets"
# "mecc"
# "newfoundland_omniquad"
# "omni_comb"
# "omniexpress_exomechip"
# "oncoarray_to_usc"
# "plco_3"
# "reach"
# "ukbiobank"  


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
"ukbiobank"
do
    for chr in {1..22}
        do echo ${OUT}/tmp/${batch}_backbone_chr${chr}
    done > ${OUT}/mergelist_$batch.txt

    merge_backbone $batch

done

