#!/bin/bash
#SBATCH -J gem_smcabc_at_three_conditions               # sensible name for the job
#SBATCH --nodes=1                    
#SBATCH -c 20
#SBATCH -t 100:00:00             # Upper time limit for the job
#SBATCH --array=0-7
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jakob.p.pettersen@ntnu.no
#SBATCH -p CPUQ
#SBATCH --account=share-nv-ibt

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
module purge
module load Anaconda3/2020.07
source /cluster/apps/eb/software/Anaconda3/2020.07/etc/profile.d/conda.sh
conda activate etcFBA
python gem_smcabc_at_three_conditions_run.py &> "../results/permuted_smcabc_res/gem_smcabc_at_three_conditions_$SLURM_ARRAY_TASK_ID.log"