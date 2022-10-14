#!/bin/bash
#SBATCH -J gem_fva_at_two_conditions               # sensible name for the job
#SBATCH --output=gem_fva_at_two_conditions.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH -t 10:00:00             # Upper time limit for the job
#SBATCH --mem=120G
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jakob.p.pettersen@ntnu.no
#SBATCH -p CPUQ
#SBATCH --account=nv-ibt
#SBATCH --export=NONE

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
module purge
module load Anaconda3/2020.07
source ~/.bash_profile
conda activate etcFBA
python gem_fva_at_two_conditions.py &> "../results/reduced_smcabc_res/gem_fva_at_two_conditions.log"
