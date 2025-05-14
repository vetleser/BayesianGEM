#!/bin/bash
#SBATCH --job-name="gem_smcsa"               # sensible name for the job
#SBATCH --output=gem_smcsa.out
#SBATCH --nodes=1                    
#SBATCH -c 20
#SBATCH -t 150:00:00             # Upper time limit for the job
#SBATCH --array=0-3
#SBATCH --mem=100G
#SBATCH --export=NONE
##SBATCH --nodelist=pompey,crassus,antony,caesar

echo "Running SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID on $(hostname)"


WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 gem_smcsa_run.py &> "../results/sa/gem_smcsa_$SLURM_ARRAY_TASK_ID.log"
