#!/bin/bash
#SBATCH -J evo_tsne               # sensible name for the job
#SBATCH --output=evo_tsne.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH --mem=60G
#SBATCH -t 10:00:00             # Upper time limit for the job
#SBATCH --array=0-9
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
python evo_tsne.py &> "../results/evo_tsne_res/evo_tsne_$SLURM_ARRAY_TASK_ID.log"
