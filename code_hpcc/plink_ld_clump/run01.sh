#!/bin/bash

#sbatch --output=corect_oncoarray_convert_plink.log --export=batch='corect_oncoarray' j01.sh
#sbatch --export=batch='corect_oncoarray' j01.sh

sbatch --export=batch='corect_oncoarray' corect_oncoarray_controls.sh

