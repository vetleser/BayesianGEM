#!/bin/bash
#SBATCH --job-name="record_flux_analyze"               # sensible name for the job
#SBATCH --output=record_flux_analyze.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH -t 00:10:00             # Upper time limit for the job
#SBATCH --mem=100G


WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}

/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 record_flux_analyze.py &> "../results/analysis/record_flux_analyze.log"