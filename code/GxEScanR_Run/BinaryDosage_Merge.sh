#!/bin/bash
# run GxEScanR

OUT=/staging/dvc/andreeki/figi
GIT=/auto/pmd-01/andreeki/FIGI_GxE_FullRun
REF=/auto/pmd-02/figi/HRC_VCF
REFBD=/auto/pmd-02/figi/HRC_BDose

cd ${OUT}

bdose_merge () {
echo "#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --job-name=bdose_merge_$1
#SBATCH --output=./logs/bdose_merge_$1.log

Rscript BinaryDosage_Merge.R $1

" | sbatch
}

# Run function here
bdose_merge 20