#!/bin/bash

for i in {1..22}
do
    awk '{print $1 "\t" $2":"$6":"$5 "\t" $3 "\t" $4 "\t" $5 "\t" $6 }' corect_oncoarray_chr${i}.bim > corect_oncoarray_chr${i}_fix.bim

done