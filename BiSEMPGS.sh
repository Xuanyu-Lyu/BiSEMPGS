#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=80gb
#SBATCH --array=1-50%50
#SBATCH -J BiSEMPGS
#SBATCH --chdir /projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS
#SBATCH --exclude bmem-rico1
#SBATCH -o %x.out%A
#SBATCH -e %x.err%A

source /curc/sw/anaconda3/latest
conda activate /projects/lessem/software/anaconda/envs/R-latest

# Index to run different parts of the script, e.g., processing different chromosomes; Here is the different simulation conditions I want to run
SIM=${SLURM_ARRAY_TASK_ID}

Rscript rc_run2simulate.R ${SIM} 
