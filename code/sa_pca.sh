#!/bin/bash
#SBATCH -J sa_pca               # sensible name for the job
#SBATCH --output=sa_pca.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -t 01:00:00             # Upper time limit for the job
#SBATCH --export=NONE

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 sa_pca.py &> "../results/sa/sa_pca.log"
