#!/bin/bash

#sbatch --output=./logs/j01_reach_chr%a.log --export=batch=reach j01.sh

sbatch --output=./logs/j01_UKB2_chr%a.log --export=batch=UKB2 j01.sh
