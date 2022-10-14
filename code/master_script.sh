#!/bin/bash

# This scripts tried to emulate the process of running all computational scripts in
# order, but is too slow too be used in practice

conda activate etcFBA

python gem_smcabc_at_three_conditions_prepare.py

for SLURM_ARRAY_TASK_ID in {0..7}
do
    python gem_smcabc_at_three_conditions_run.py &> "../results/permuted_smcabc_res/gem_smcabc_at_three_conditions_$SLURM_ARRAY_TASK_ID.log"
done

python gem_smcabc_tournament_prepare.py


for SLURM_ARRAY_TASK_ID in {0..7}
do 
    python gem_smcabc_tournament.py &> "../results/evo_tournament/gem_smcevo_tournament_$SLURM_ARRAY_TASK_ID.log"
done

python gem_smcabc_truncation_prepare.py

for SLURM_ARRAY_TASK_ID in {0..7}
do 
    python gem_smcabc_truncation.py &> "../results/evo_truncation/gem_smcevo_truncation_$SLURM_ARRAY_TASK_ID.log"
done


python sample_pca.py &> "../results/permuted_smcabc_res/sample_pca.log"

python evo_pca.py &> "../results/evo_pca.log"

python reduced_pca.py &> "../results/reduced_smcabc_res/reduced_pca.log"


python gem_fva_at_three_conditions.py &> "../results/permuted_smcabc_res/gem_fva_at_three_conditions.log"

python gem_fva_at_two_conditions.py &> "../results/reduced_smcabc_res/gem_fva_at_two_conditions.log"

python evo_fva.py &> "../results/evo_fva.log"

# python code/compute_particle_distances.py

python reduce_data_size.py

echo "Work flow done"
