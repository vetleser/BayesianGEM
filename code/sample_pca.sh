#!/bin/bash
#SBATCH -J sample_pca               # sensible name for the job
#SBATCH --output=sample_pca.out
#SBATCH --nodes=1                    
#SBATCH -c 8
#SBATCH --mem=100G
#SBATCH -t 10:00:00             # Upper time limit for the job
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=vetleser@stud.ntnu.no
#SBATCH -p CPUQ
#SBATCH --account=nv-ibt
#SBATCH --export=NONE

WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}
module purge
module load Anaconda3/2020.07
source ~/.bash_profile
conda activate etcFBA
python sample_pca.py &> "../results/permuted_smcabc_res/sample_pca.log"
