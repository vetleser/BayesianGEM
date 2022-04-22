#!/bin/bash
#SBATCH -J gem_smcevo_at_three_conditions               # sensible name for the job
#SBATCH --output=gem_smcevo_at_three_conditions.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH -t 00:01:00             # Upper time limit for the job
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jakob.p.pettersen@ntnu.no
#SBATCH -p CPUQ
#SBATCH --account=nv-ibt

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
module purge
module load Anaconda3/2020.07
source /cluster/apps/eb/software/Anaconda3/2020.07/etc/profile.d/conda.sh
conda activate etcFBA
python gem_smcevo_at_three_conditions.py &> "../results/gem_smcevo_at_three_conditions.log"
