#!/bin/bash
#SBATCH --mem=6G
#SBATCH -n 1
#SBATCH -N 1

ml R

echo "Starting R script..."

# [ASE Based Data] [General Data] [ASE File Start] [ASE File stop] [Gen File Start] [Gen File Stop] [Iteration]

Rscript randGeneCorrelationLinearModel_v_1.7.r $1 $2 $3 $4 $5 $6 $7
echo "Complete!"
