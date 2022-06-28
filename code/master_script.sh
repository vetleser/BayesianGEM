#!/bin/bash

# This scripts tried to emulate the process of running all computational scripts in
# order, but is too slow too be used in practice

conda activate etcFBA

python gem_smcabc_at_three_conditions_prepare.py

for SLURM_ARRAY_TASK_ID in {0..7}
do
    python gem_smcabc_at_three_conditions_run.py &> "../results/permuted_smcabc_res/gem_smcabc_at_three_conditions_$SLURM_ARRAY_TASK_ID.log"
done

python gem_smcevo_at_three_conditions.py &> "../results/gem_smcevo_at_three_conditions.log"

python sample_pca.py &> "../results/permuted_smcabc_res/sample_pca.log"

python evo_pca.py &> "../results/evo_pca.log"

python reduced_pca.py &> "../results/permuted_smcabc_res/reduced_pca.log"

python evo_tsne_prepare.py

for SLURM_ARRAY_TASK_ID in {0..9}
do
    python evo_tsne.py &> "../results/evo_tsne_res/evo_tsne_$SLURM_ARRAY_TASK_ID.log"
done

python gem_fva_at_three_conditions.py &> "../results/permuted_smcabc_res/gem_fva_at_three_conditions.log"

python evo_fva.py &> "../results/evo_fva.log"

python code/compute_particle_distances.py

python reduce_data_size.py

echo "Work flow done"