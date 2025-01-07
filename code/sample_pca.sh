#!/bin/bash
#SBATCH -J sample_pca               # sensible name for the job
#SBATCH --output=sample_pca.out
#SBATCH --nodes=1                    
#SBATCH -c 8
#SBATCH --mem=100G
#SBATCH -t 10:00:00             # Upper time limit for the job
#SBATCH --export=NONE

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 sample_pca.py &> "../results/permuted_smcabc_res/sample_pca.log"
