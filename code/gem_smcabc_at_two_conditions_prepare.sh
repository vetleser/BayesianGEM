#!/bin/bash
#SBATCH --job-name="gem_smcabc_at_two_conditions_prepare"               # sensible name for the job
#SBATCH --output=gegem_smcabc_at_two_conditions_prepare.out
#SBATCH --nodes=1                    
#SBATCH -c 20
#SBATCH -t 10:00:00             # Upper time limit for the job
#SBATCH --mem=100G
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=vetleser@stud.ntnu.no
#SBATCH --export=NONE

/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 gem_smcabc_at_two_conditions_prepare.py