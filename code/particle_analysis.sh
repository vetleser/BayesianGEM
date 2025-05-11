#!/bin/bash
#SBATCH -J particle_analysis               # sensible name for the job
#SBATCH --output=particle_analysis.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -t 01:00:00             # Upper time limit for the job
#SBATCH --export=NONE

export HOME=/triumvirate/home/vetleser  # Set the HOME environment variable explicitly


WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 particle_analysis.py &> "../results/analysis/particle_analysis.log"
