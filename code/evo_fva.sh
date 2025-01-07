#!/bin/bash
#SBATCH -J evo_fva               # sensible name for the job
#SBATCH --output=evo_fva.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH -t 10:00:00             # Upper time limit for the job
#SBATCH --mem=120G
#SBATCH --export=NONE

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
/triumvirate/home/vetleser/.conda/envs/etcFBA_v2/bin/python3 evo_fva.py &> "../results/crowdingDE/evo_fva.log"
