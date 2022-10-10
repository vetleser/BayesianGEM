#!/bin/bash
#SBATCH -J evo_fva               # sensible name for the job
#SBATCH --output=evo_fva.out
#SBATCH --nodes=1                    
#SBATCH -c 20
#SBATCH -t 10:00:00             # Upper time limit for the job
#SBATCH --mem=120G
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
python evo_fva.py &> "../results/evo_fva.log"
