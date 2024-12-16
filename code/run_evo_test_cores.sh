#!/bin/bash
#SBATCH --job-name="run_evo_test_cores"               # sensible name for the job
#SBATCH --output=run_evo_test_cores.out
#SBATCH --nodes=1
#SBATCH -c 35 #ENDRE HER?
#SBATCH -t 05:00:00             # Upper time limit for the job
#SBATCH --array=0-1
#SBATCH --mem=100G
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=vetleser@stud.ntnu.no
#SBATCH --export=NONE



WORKDIR=${SLURM_SUBMIT_DIR}
cd ${WORKDIR}

/triumvirate/home/vetleser/.conda/envs/etcFBA/bin/python3 run_evo_test_cores.py &> "../results/test/run_evo_test_cores_$SLURM_ARRAY_TASK_ID.log"