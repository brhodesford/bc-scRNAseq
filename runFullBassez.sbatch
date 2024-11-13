#!/bin/bash
#SBATCH -J integration
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256g
#SBATCH -t 0-05:00:00
#SBATCH --mail-user=brford
#SBATCH --mail-type=END

# Load modules
module load CBI r

#Change directory
cd /c4/home/brford/breastCancerProject/Bassez/

# Run R script with provided file and arguments
Rscript runFullBassezData.R
