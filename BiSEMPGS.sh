#!/bin/bash
#SBATCH --qos=blanca -preemptable
#SBATCH --time=12:00:00
#SBATCH --ntasks=2
#SBATCH --mem=16gb
#SBATCH --array=1-10
#SBATCH -J BiSEMPGS
#SBATCH --chdir /projects/xuly4739/BiSEMPGS/BiSEMPGS
#SBATCH --exclude bmem-rico1,bnode010[1-5]
#SBATCH -o %x.out%A
#SBATCH -e %x.err%A
source /curc/sw/anaconda3/R-4.3.3
conda activate /projects/lessem/software/anaconda/envs/R-4.3.3

# Index to run different parts of the script, e.g., processing different chromosomes; Here is the different simulation conditions I want to run
SIM=${SLURM_ARRAY_TASK_ID}

Rscript rc_run2simulate.R ${CHR} 