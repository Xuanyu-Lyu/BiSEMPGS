#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --time=20:00:00
#SBATCH --ntasks=16
#SBATCH --mem=40gb
#SBATCH -J BiSEMPGS_fit_lowtol
#SBATCH --chdir /projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS
#SBATCH --exclude bmem-rico1
#SBATCH -o %x.out%A
#SBATCH -e %x.err%A

source /curc/sw/anaconda3/latest
conda activate /projects/lessem/software/anaconda/envs/R-4.4.0

# Index to run different parts of the script, e.g., processing different chromosomes; Here is the different simulation conditions I want to run
#SIM=${SLURM_ARRAY_TASK_ID}

Rscript BiSEMPGS_fit_rc.R 