#!/bin/bash

# SLURM job scheduling directives for HPC cluster resource allocation
#SBATCH --account=def-bureau
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-01:30


# Output and notification configuration
#SBATCH --output=/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/paper3_Analyzes/Merlin_Recomb/Slurm_%A.out
#SBATCH --mail-user=samir.oubninte.1@ulaval.ca
#SBATCH --mail-type=FAIL       # Send email notification only on job failure


# This script performs Merlin-based recombination rate analysis .
# It is designed to run on an HPC cluster using SLURM job scheduler and executes R-based analysis



# Example of how to submit this script with chromosome parameter via loop:
# sbatch for chrom in 22; do sbatch --export=ALL,chrom=${chrom} Data_processing.sh ; done

# Display status message to indicate module loading phase
echo "**************************** Import module ****************************"

# Load required software modules for the HPC environment
module load StdEnv/2020
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14
module load r/4.2.1

# Execute R scripts for Merlin recombination analysis
#Infer recombination events for non‑split families based on Merlin output
Rscript /home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/paper3_Analyzes/Recomb_Merlin_avantDecoup.R  
#Infer recombination events for splited families based on Merlin output
Rscript /home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/paper3_Analyzes/Recomb_Merlin.R 

