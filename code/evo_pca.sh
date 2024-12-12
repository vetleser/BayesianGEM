#!/bin/bash
#SBATCH -J evo_pca               # sensible name for the job
#SBATCH --output=evo_pca.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -t 01:00:00             # Upper time limit for the job
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=vetleser@stud.ntnu.no
#SBATCH -p CPUQ
#SBATCH --account=share-nv-ibt
#SBATCH --export=NONE

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
module purge
module load Anaconda3/2020.07
source ~/.bash_profile
conda activate etcFBA
python evo_pca.py &> "../results/crowdingDE/evo_pca.log"
