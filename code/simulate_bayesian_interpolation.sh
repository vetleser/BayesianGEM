#!/bin/bash
#SBATCH -J simulate_bayesian_interpolation               # sensible name for the job
#SBATCH --output=simulate_bayesian_interpolation.out
#SBATCH --nodes=1                    
#SBATCH -c 20
#SBATCH -t 100:00:00             # Upper time limit for the job
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
python simulate_bayesian_interpolation.py &> "../results/permuted_smcabc_res/simulate_bayesian_interpolation.log"
