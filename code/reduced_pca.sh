#!/bin/bash
#SBATCH -J reduced_pca               # sensible name for the job
#SBATCH --output=reduced_pca.out
#SBATCH --nodes=1                    
#SBATCH -c 4
#SBATCH --mem=100G
#SBATCH -t 10:00:00             # Upper time limit for the job
#SBATCH --export=NONE

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 reduced_pca.py &> "../results/reduced_smcabc_res/reduced_pca.log"
