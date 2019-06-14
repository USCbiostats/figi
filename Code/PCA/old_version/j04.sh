#!/bin/bash
# extract same backbone SNPs from Thousand Genome files

REF=/auto/pmd-02/figi/andreeki/Refs/1KGP
OUT=/staging/dvc/andreeki/PCA

cd ${OUT}

merge_backbone_kgp () {

echo "#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --job-name=j04_kgp_$1
#SBATCH --output=./logs/j04_kgp_$1.log

cd ${OUT}

plink --bfile ${REF}/kgp.chr$1.biallelic --extract ${OUT}/FIGI_PC_Backbone_Sample_30K_rsID.txt --memory 16000 --make-bed --out ${OUT}/tmp/kgp_backbone_chr$1

" | sbatch
}

for chr in {1..22}
do
    merge_backbone_kgp $chr
done



