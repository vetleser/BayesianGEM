#!/bin/bash
#SBATCH --job-name="gem_smcabc_at_three_conditions_resume"               # sensible name for the job
#SBATCH --output=gem_smcabc_at_three_conditions_resume.out
#SBATCH --nodes=1                    
#SBATCH -c 20
#SBATCH -t 00:30:00             # Upper time limit for the job
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
python gem_smcabc_at_three_conditions_resume.py &> "../results/permuted_smcabc_res/gem_smcabc_at_three_conditions_resume.log"
