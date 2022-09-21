#!/bin/bash
#SBATCH --job-name="slurm_debug"               # sensible name for the job
#SBATCH --output=slurm_debug.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH -t 00:01:00             # Upper time limit for the job
#SBATCH --array=0-15
#SBATCH --mem=50G
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
python slurm_debug.py &> "slurm_debug_$SLURM_ARRAY_TASK_ID.log"

