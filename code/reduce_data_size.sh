#!/bin/bash
#SBATCH -J reduce_data_size               # sensible name for the job
#SBATCH --output=reduce_data_size.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -t 01:00:00             # Upper time limit for the job
#SBATCH --export=NONE

export HOME=/triumvirate/home/vetleser  # Set the HOME environment variable explicitly


WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
/triumvirate/home/vetleser/.conda/envs/etcFBA_v2/bin/python3 reduce_data_size.py &> "../results/crowdingDE/reduce_data_size.log"
