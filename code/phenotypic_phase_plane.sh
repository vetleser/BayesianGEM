#!/bin/bash
#SBATCH -J phenotypic_phase_plane               # sensible name for the job
#SBATCH --output=phenotypic_phase_plane.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH -t 10:00:00             # Upper time limit for the job
#SBATCH --mem=120G


WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 phenotypic_phase_plane.py &> "../results/analysis/phenotypic_phase_plane.log"