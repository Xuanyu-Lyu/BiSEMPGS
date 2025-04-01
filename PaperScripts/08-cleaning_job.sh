#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --time=6:00:00
#SBATCH --ntasks=4
#SBATCH --mem=40gb
#SBATCH -J BiSEMPGS_dataprep
#SBATCH --array=1-100%35
#SBATCH --chdir /projects/xuly4739/R-Projects/BiSEMPGS/BiSEMPGS
#SBATCH --exclude bmem-rico1
#SBATCH -o %x.out%A
#SBATCH -e %x.err%A

source /curc/sw/anaconda3/latest
conda activate /projects/lessem/software/anaconda/envs/R-latest

# Index to run different parts of the script, e.g., processing different chromosomes; Here is the different simulation conditions I want to run
SIM=${SLURM_ARRAY_TASK_ID}

Rscript PaperScripts/08-dataPrep.R
#Rscript PaperScripts/03-createMoreTxt.R
