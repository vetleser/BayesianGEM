#!/bin/bash
#SBATCH --job-name="gem_smcabc_at_two_conditions"               # sensible name for the job
#SBATCH --output=gem_smcabc_at_two_conditions.out
#SBATCH --nodes=1                    
#SBATCH -c 20
#SBATCH -t 100:00:00             # Upper time limit for the job
#SBATCH --array=0-3
#SBATCH --mem=50G
#SBATCH --export=NONE

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 gem_smcabc_at_two_conditions.py &> "../results/reduced_smcabc_res/gem_smcabc_at_two_conditions_$SLURM_ARRAY_TASK_ID.log"
