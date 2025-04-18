#!/bin/bash
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=300gb
#SBATCH --array=1-100%20
#SBATCH -J BiSEMPGS_simulation_new
#SBATCH --chdir /projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS
#SBATCH -o %x.out%A
#SBATCH -e %x.err%A

source /curc/sw/anaconda3/latest
conda activate /projects/lessem/software/anaconda/envs/R-latest

# Index to run different parts of the script, e.g., processing different chromosomes; Here is the different simulation conditions I want to run
SIM=${SLURM_ARRAY_TASK_ID}

Rscript PaperScripts/01-run2simulate.R ${SIM} 
