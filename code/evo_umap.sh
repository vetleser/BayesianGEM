#!/bin/bash
#SBATCH -J evo_umap              # sensible name for the job
#SBATCH --output=evo_umap.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH -t 01:00:00             # Upper time limit for the job
#SBATCH --mem=120G
#SBATCH --export=NONE

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 evo_umap.py &> "../results/crowdingDE/evo_umap.log"