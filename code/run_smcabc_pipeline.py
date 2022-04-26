#!/usr/bin/env python
# coding: utf-8
# Master SLURM script to control job dependencies for smcabc permutation pipeline
import subprocess

# submit the first job
cmd = "sbatch gem_smcabc_at_three_conditions_run.sh"
print(f"Submitting fitting job with command: {cmd}")
status, main_jobnum = subprocess.getstatusoutput(cmd)
if (status == 0 ):
    print(f"Fitting job is {main_jobnum}")
else:
    print("Error submitting fitting job")

dependent_jobs = {'PCA job': 'sample_pca.sh', 'FVA job': 'gem_fva_at_three_conditions.sh', 'Interpolation job': 'simulate_bayesian_interpolation.sh'}
for jobname, script in dependent_jobs.items():
    cmd = f"sbatch --depend=afterok:{main_jobnum} {script}"
    status,jobnum = subprocess.getstatusoutput(cmd)
    if (status == 0 ):
        print(f"{jobname} is {jobnum}")
    else:
        print(f"Error submitting {jobname}")
