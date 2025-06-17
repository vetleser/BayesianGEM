#!/bin/bash
#SBATCH -J create_dataframe_sa               # sensible name for the job
#SBATCH --output=create_dataframe_sa.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -t 01:00:00             # Upper time limit for the job
#SBATCH --export=NONE

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 create_dataframe_sa.py &> "../results/sa/create_dataframe_sa.log"
