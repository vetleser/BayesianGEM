#!/bin/bash
#SBATCH --job-name="test_record_flux"               # sensible name for the job
#SBATCH --output=test_record_flux.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH -t 00:10:00             # Upper time limit for the job
#SBATCH --mem=100G


WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}

/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 test_record_flux.py &> "../results/analysis/test_record_flux.log"