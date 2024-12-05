#!/bin/bash
#SBATCH --job-name="gem_smcevo"               # sensible name for the job
#SBATCH --output=gem_smcevo.out
#SBATCH --nodes=1                    
#SBATCH -c 20
#SBATCH -t 50:00:00             # Upper time limit for the job
#SBATCH --array=0-11
#SBATCH --mem=100G
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=vetleser@stud.ntnu.no
#SBATCH --export=NONE

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}

/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 gem_smcevo_run.py &> "../results/crowdingDE/gem_smcevo_$SLURM_ARRAY_TASK_ID.log"
