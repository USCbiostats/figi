#!/bin/bash

# ----------------------------------------------------------------------------#
# Preparing chr X imputation for victor's data
#
# Before imputation of chr X, need to:
# 1) check-sex
# 2) fix heterozygous male coding
#
# important - already ran wrayner script on victor's data
# use those as a starting point (FIGI_tmp-updated-chr23)
#
#-----------------------------------------------------------------------------#


# 1) sex-check
#plink --bfile FIGI_tmp-updated-chr23 --check-sex --out vm_tmp1
plink --bfile FIGI_tmp-updated-chr23 --indep-pairphase 20000 2000 0.5 --make-bed --out TMP1
plink --bfile FIGI_tmp-updated-chr23 --exclude TMP1.prune.out --make-bed --out TMP2
plink --bfile FIGI_tmp-updated-chr23 --check-sex --out TMP2

# manually fix samples
# switch sexes for the following:
#   VM1261                              S16759            2            1      PROBLEM       0.9013
#   VM1330                              S16828            2            1      PROBLEM       0.9996
#   VM1331                              Q16829            1            2      PROBLEM     -0.09125
#   VM1374                              J16872            1            0      PROBLEM       0.4226

plink --bfile FIGI_tmp-updated-chr23 --update-sex vm_sexchange.txt --make-bed --out TMP3
plink --bfile TMP3 --set-hh-missing --make-bed --out FIGI_tmp-updated-chr23_HRC 
plink --bfile FIGI_tmp-updated-chr23_HRC  --recode vcf --out FIGI_tmp-updated-chr23_HRC 

awk -F'\t' -v OFS='\t' '{if ($1 ~ /^[^#]/) {gsub(/23/,"X",$1); print} else {print}}' FIGI_tmp-updated-chr23_HRC.vcf > FIGI_tmp-updated-chr23_HRC-fix.vcf
vcf-sort FIGI_tmp-updated-chr23_HRC-fix.vcf | bgzip -c > FIGI_tmp-updated-chr23_HRC-fix.vcf.gz

