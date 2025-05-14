#!/bin/bash
#SBATCH --job-name="gem_smcsa_prepare"               # sensible name for the job
#SBATCH --output=gem_smcsa_prepare.out
#SBATCH --nodes=1                    
#SBATCH -c 20
#SBATCH -t 10:00:00             # Upper time limit for the job
#SBATCH --mem=100G

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}


/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 gem_smcsa_prepare.py &> "../results/sa/gem_smcsa_prepare.log"