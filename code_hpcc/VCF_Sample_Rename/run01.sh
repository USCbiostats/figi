#!/bin/bash

sbatch --output=./logs/VCF_sample_rename_UKB2_chr%a.log --export=batch=UKB2 job01.sh
