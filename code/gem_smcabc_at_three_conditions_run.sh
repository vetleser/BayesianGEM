#!/bin/bash
#SBATCH -J gem_smcabc_at_three_conditions               # sensible name for the job
#SBATCH -N 8                    
#SBATCH --ntasks-per-node=1     # 1 task per node
#SBATCH -c 20
#SBATCH -t 00:05:00             # Upper time limit for the job
#SBATCH -a 0-7
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jakob.p.pettersen@ntnu.no
#SBATCH -p CPUQ

# module purge
# module load Anaconda3/2020.07
conda activate etcFBA
python gem_smcabc_at_three_conditions_run.py &> "gem_smcabc_at_three_conditions_$SLURM_ARRAY_JOB_ID.log"
