#!/bin/bash
#SBATCH -J visualization_umap              # sensible name for the job
#SBATCH --output=visualization_umap.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH -t 01:00:00             # Upper time limit for the job
#SBATCH --mem=120G
#SBATCH --export=NONE

export HOME=/triumvirate/home/vetleser  # Set the HOME environment variable explicitly


WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
/triumvirate/home/vetleser/.conda/envs/umap/bin/python3 visualization_umap.py &> "../results/crowdingDE/visualization_umap.log"