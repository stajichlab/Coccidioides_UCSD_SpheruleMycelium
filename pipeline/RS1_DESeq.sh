#!/bin/bash
#SBATCH --nodes 1 --ntasks 1 --out logs/DEseq_RS1.log --mem 4G 
#SBATCH --time 2:00:00 -p short
module load R
module load R/3.6.0
which Rscript
mkdir -p plots reports
Rscript Rscripts/kallisto_profile_RS1.R

