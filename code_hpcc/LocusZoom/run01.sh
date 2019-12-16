#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=19-20

# LocusZoom local includes 1000G data used to calculate LD structure
# can also provide your own plink dataset, which I'm creating here based on the axiom_mecc_cfr_ky imputation batch 
# (uses hard calls, which HRC provides in the VCF file)

plink --vcf /auto/pmd-02/figi/HRC_extract/axiom_mecc_cfr_ky/chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz --keep-allele-order --make-bed --out chr${SLURM_ARRAY_TASK_ID}

cp chr${SLURM_ARRAY_TASK_ID}.bim chr${SLURM_ARRAY_TASK_ID}.bim.bck
awk '{$2="chr"$2}1' chr${SLURM_ARRAY_TASK_ID}.bim.bck > chr${SLURM_ARRAY_TASK_ID}.bim
