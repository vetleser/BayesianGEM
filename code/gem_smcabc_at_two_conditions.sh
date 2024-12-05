#!/bin/bash
#SBATCH --job-name="gem_smcabc_at_two_conditions"               # sensible name for the job
#SBATCH --output=gem_smcabc_at_two_conditions.out
#SBATCH --nodes=1                    
#SBATCH -c 20
#SBATCH -t 100:00:00             # Upper time limit for the job
#SBATCH --array=0-3
#SBATCH --mem=50G
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=vetleser@stud.ntnu.no
#SBATCH -p CPUQ
#SBATCH --account=nv-ibt
#SBATCH --export=NONE

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
module purge
module load Anaconda3/2020.07
source ~/.bash_profile
conda activate etcFBA
python gem_smcabc_at_two_conditions.py &> "../results/reduced_smcabc_res/gem_smcabc_at_two_conditions_$SLURM_ARRAY_TASK_ID.log"
