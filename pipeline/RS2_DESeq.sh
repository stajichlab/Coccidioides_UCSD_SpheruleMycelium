#!/bin/bash
#SBATCH --nodes 1 --ntasks 1 --out logs/DEseq.log --mem 4G 
#SBATCH --time 2:00:00 -p short

mkdir -p plots reports
Rscript Rscripts/kallisto_profile_RS2.R
