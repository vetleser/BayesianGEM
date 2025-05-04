#!/bin/bash
#SBATCH -J filtering_data               # sensible name for the job
#SBATCH --output=filtering_data.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -t 01:00:00             # Upper time limit for the job
#SBATCH --export=NONE

export HOME=/triumvirate/home/vetleser  # Set the HOME environment variable explicitly


WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 filtering_data.py &> "../results/analysis/filtering_data.log"
