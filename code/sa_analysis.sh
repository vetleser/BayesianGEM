#!/bin/bash
#SBATCH --job-name="sa_analysis"               # sensible name for the job
#SBATCH --output=sa_analysis.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH -t 00:10:00             # Upper time limit for the job
#SBATCH --mem=100G


WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}

/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 sa_analysis.py &> "../results/analysis/sa_analysis.log"