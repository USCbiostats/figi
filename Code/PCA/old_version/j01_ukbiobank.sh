#!/bin/bash
# Extract snp backbone (random sample) from VCFs (all batches + chromosomes)
# (This part doesn't need to be ran over and over, just once to extract VCFs)

REF=/auto/pmd-02/figi/HRC_VCF_SampleRename
OUT=/staging/dvc/andreeki/pca_ibd

cd ${OUT}

extract_backbone () {

echo "#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --job-name=j01b_ukb_$2
#SBATCH --output=./logs/j01b_ukb_$2.log

plink2 --vcf ${OUT}/tmp/$1_backbone_chr$2.recode.vcf dosage=DS --memory 16000 --double-id --keep-allele-order --make-bed --out ${OUT}/tmp/$1_backbone_chr$2

" | sbatch
}

for batch in "ukbiobank"
do
    for chr in {6..6}
    do
        extract_backbone $batch $chr
    done
done




plink2 --vcf ./tmp/ukbiobank_backbone_chr6.recode.vcf dosage=DS --memory 16000 --double-id --keep-allele-order --make-bed --out ./tmp/ukbiobank_backbone_chr6

# might need to add hard call threshold because make-bed was setting uncertain dosages to missing! Fuuuu
plink2 --vcf check.recode.vcf dosage=DS --double-id --keep-allele-order --hard-call-threshold 0.49 --make-bed --out wtf