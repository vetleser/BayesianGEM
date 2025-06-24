#!/bin/bash
#SBATCH --job-name="analyze_log_files"               # sensible name for the job
#SBATCH --output=analyze_log_files.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH -t 05:00:00             # Upper time limit for the job
#SBATCH --mem=100G


WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}

/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 analyze_log_files.py &> "../results/analysis/analyze_log_files.log"