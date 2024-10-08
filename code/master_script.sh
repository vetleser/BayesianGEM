#!/bin/bash

# This scripts tried to emulate the process of running all computational scripts in
# order, but is too slow too be used in practice

conda activate etcFBA

python gem_smcabc_at_three_conditions_prepare.py

for SLURM_ARRAY_TASK_ID in {0..7}
do
    python gem_smcabc_at_three_conditions_run.py &> "../results/permuted_smcabc_res/gem_smcabc_at_three_conditions_$SLURM_ARRAY_TASK_ID.log"
done

python gem_smcevo_prepare.py


for SLURM_ARRAY_TASK_ID in {0..7}
do 
    python gem_smcevo_run.py &> "../results/CrowdingDE/gem_smcevo_$SLURM_ARRAY_TASK_ID.log"
done

python toy_example.py

python sample_pca.py &> "../results/permuted_smcabc_res/sample_pca.log"

python evo_pca.py &> "../results/crowdingDE/evo_pca.log"

python reduced_pca.py &> "../results/reduced_smcabc_res/reduced_pca.log"


python gem_fva_at_three_conditions.py &> "../results/permuted_smcabc_res/gem_fva_at_three_conditions.log"

python gem_fva_at_two_conditions.py &> "../results/reduced_smcabc_res/gem_fva_at_two_conditions.log"

python evo_fva.py &> "../results/crowdingDE/evo_fva.log"

python reduce_data_size.py

echo "Workflow done"
