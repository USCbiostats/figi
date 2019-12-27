#!/bin/bash

sbatch --job-name=clump --output=/staging/dvc/andreeki/clump_results/logs/plink_asp_ref_ldclump_chiSqGxE_chr%a.log --export=exposure=asp_ref plink_ld_clump.sh


