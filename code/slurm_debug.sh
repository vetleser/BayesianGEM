#!/bin/bash
#SBATCH --job-name="slurm_debug"               # sensible name for the job
#SBATCH --output=slurm_debug.out
#SBATCH --nodes=1                    
#SBATCH -c 1
#SBATCH -t 50:00:00             # Upper time limit for the job
#SBATCH --array=0-15
#SBATCH --mem=50G
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jakob.p.pettersen@ntnu.no
#SBATCH -p CPUQ
#SBATCH --account=share-nv-ibt
#SBATCH --export=NONE
