#!/bin/bash
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH -N 1

ml R

echo "Starting R script..."

# [ASE Based Data] [General Data] [ASE File Start] [ASE File stop] [Gen File Start] [Gen File Stop]

Rscript GeneCorrelationLinearModel_v_1.6.r $1 $2 $3 $4 $5 $6
echo "Complete!"
