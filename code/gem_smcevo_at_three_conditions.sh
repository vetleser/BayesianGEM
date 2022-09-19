#!/bin/bash
#SBATCH --job-name="gem_smcevo_at_three_conditions"               # sensible name for the job
#SBATCH --output=gem_smcevo_at_three_conditions.out
#SBATCH --nodes=1                    
#SBATCH -c 20
#SBATCH -t 100:00:00             # Upper time limit for the job
#SBATCH --array=0-7
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jakob.p.pettersen@ntnu.no
#SBATCH -p CPUQ
#SBATCH --account=share-nv-ibt
#SBATCH --export=NONE

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
module purge
module load Anaconda3/2020.07
source ~/.bash_profile
conda activate etcFBA
python gem_smcevo_at_three_conditions_run.py &> "../results/permuted_smcevo_res/gem_smcevo_at_three_conditions_$SLURM_ARRAY_TASK_ID.log"
