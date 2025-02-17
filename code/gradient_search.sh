#!/bin/bash
#SBATCH -J gradient_search               # sensible name for the job
#SBATCH --output=gradient_search.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --array=0-9
#SBATCH -t 01:00:00             # Upper time limit for the job
#SBATCH --export=NONE

export HOME=/triumvirate/home/vetleser  # Set the HOME environment variable explicitly


WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 gradient_search.py &> "../results/analysis/gradient_search_$SLURM_ARRAY_TASK_ID.log"
