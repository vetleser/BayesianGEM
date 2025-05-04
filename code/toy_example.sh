#!/bin/bash
#SBATCH --job-name="toy_example"               # sensible name for the job
#SBATCH --output=toy_example.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH -t 00:10:00             # Upper time limit for the job
#SBATCH --mem=100G


WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}

/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 toy_example.py &> "../results/sa/toy_example_2.log"