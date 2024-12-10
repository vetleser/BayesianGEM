#!/bin/bash
#SBATCH --job-name="test_gem_smcevo"               # sensible name for the job
#SBATCH --output=test_gem_smcevo.out
#SBATCH --nodes=1                    
#SBATCH -c 20
#SBATCH -t 01:00:00             # Upper time limit for the job
#SBATCH --array=0-11
#SBATCH --mem=100G
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=vetleser@stud.ntnu.no
#SBATCH --export=NONE

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}

/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 test_evo_etc.py &> "../results/test/gem_smcevo_$SLURM_ARRAY_TASK_ID.log"