#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti

# repeat this step just to be thorough - get sample names from ALL vcf files
# to ensure that sample names match all of the sample names in samplefile from Yi (v2.1)

for chr in {1..22}
do
    bcftools query -l /auto/pmd-02/figi/HRC_VCF_SampleRename/${batch}_chr${chr}.recode.vcf > sample_list_vcf_${batch}_chr${chr}.txt