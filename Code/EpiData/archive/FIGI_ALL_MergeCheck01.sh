#!/bin/bash

# For each batch of imputation:
# - check that every chromosome has exact same sample list
# ** this is because oncoarray_to_usc did not have same samples in every chromosome (some missing one, some not)

OUT=/staging/dvc/andreeki/bdose/qc
REFBD=/auto/pmd-02/figi/HRC_BDose
GIT=/auto/pmd-01/andreeki/FIGI_EpiData

cd ${OUT}

check_bdose_files () {
echo "#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --mail-type=END
#SBATCH --job-name=bdose_qc_$1_$2
#SBATCH --output=./qc/bdose_qc_$1_$2.log

Rscript ${GIT}/FIGI_ALL_MergeCheck01.R $1 $2

" | sbatch
}

for chr in {1..22}
do
    check_bdose_files axiom_acs_aus_nf ${chr}
done

