#!/bin/bash
#SBATCH --job-name="visualization_toy_example"               # sensible name for the job
#SBATCH --output=visualization_toy_example.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH -t 01:00:00             # Upper time limit for the job
#SBATCH --mem=100G


WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}

/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 visualization_toy_example.py &> "../results/toy_example/visualization_toy_example.log"