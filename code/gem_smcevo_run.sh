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
python gem_smcevo_run.py &> "../results/crowdingDE/gem_smcevo_$SLURM_ARRAY_TASK_ID.log"
