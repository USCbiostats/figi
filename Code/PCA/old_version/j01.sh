#!/bin/bash
# Extract snp backbone (random sample) from VCFs (all batches + chromosomes)
# (This part doesn't need to be ran over and over, just once to extract VCFs)

REF=/auto/pmd-02/figi/HRC_VCF_SampleRename
OUT=/staging/dvc/andreeki/PCA

cd ${OUT}

extract_backbone () {

echo "#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --job-name=j01_$1_$2
#SBATCH --output=./logs/j01_$1_$2.log

vcftools --gzvcf ${REF}/$1_chr$2.vcf.gz --positions ${OUT}/FIGI_PC_Backbone_Sample_30K.txt --recode --out ${OUT}/tmp/$1_backbone_chr$2

plink --vcf ${OUT}/tmp/$1_backbone_chr$2.recode.vcf --double-id --snps-only --biallelic-only --keep-allele-order --make-bed --out ${OUT}/tmp/$1_backbone_chr$2

" | sbatch
}

# loop over all batches and chromosomes..

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
"axiom_mecc_cfr_ky"
do
    for chr in {1..22}
    do
        extract_backbone $batch $chr
    done
done

