#!/bin/bash

#sbatch --output=./logs/j02_reach.log --export=batch=reach j02.sh
#sbatch --output=./logs/j02_UKB1.log --export=batch=UKB1 j02.sh
sbatch --output=./logs/j02_UKB2.log --export=batch=UKB2 j02.sh
