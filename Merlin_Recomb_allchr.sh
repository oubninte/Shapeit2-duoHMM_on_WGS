#!/bin/bash

#SBATCH --account=def-bureau
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-01:30

#SBATCH --output=/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/paper3_Analyzes/Merlin_Recomb/Slurm_%A.out
#SBATCH --mail-user=samir.oubninte.1@ulaval.ca
#SBATCH --mail-type=FAIL

# to run cmd example :
  #sbatch for chrom in 22; do sbatch --export=ALL,chrom=${chrom} Data_processing.sh ; done
  
  



echo "**************************** Import module ****************************"
module load StdEnv/2020
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14
module load r/4.2.1



Rscript /home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/paper3_Analyzes/Recomb_Merlin_avantDecoup.R
Rscript /home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/paper3_Analyzes/Recomb_Merlin.R