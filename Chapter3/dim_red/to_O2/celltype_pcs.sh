#!/bin/bash

#SBATCH -c 1
#SBATCH -t 2-00:00
#SBATCH -p medium
#SBATCH --mem=64GB
#SBATCH -o celltype_pcs_%j.out
#SBATCH -e celltype_pcs_%j.err

module load gcc/6.2.0
module load R/4.1.1
export R_LIBS_USER="~/R/4.1.1/library"
srun Rscript celltype_pcs.R

