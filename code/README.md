## Description
This folder contains the scripts that carry out all analysis reported in the paper. Some scripts assumes the job is run under the SLURM workload manager. While these scripts will require some modification to run under your perferred computing environment, such modifications should be relatively feasible to carry out.

## Backend scripts
* `etcpy` - the scripts for incorporating temperature parameters into the enzyme constrained Yeast 7.6 model, the detail introduction can be found from `etcpy/README.md`.
* `GEMS.py` - the functionality for simulating the aerobic, anaerobic growth rate in batch cultivation and aerobic fluxes in chemostat cultivation, as well as a list of distance functions used for SMC-ABC and the evolutionary approach. It also includes corresponding FVA functionality for the batch and chemostat cultivations.
* `abc_etc.py` - the functionality to perform SMC-ABC approach.
* `evo_etc` - Like `abc_etc.py`, but uses an evolutionary algorithm instead of SMC-ABC
* `permute_parameters.py` - Utility script for shuffling parameters and thus creating permuted priors
* `random_sampler.py` - Utility script for creating random numbers for SMC-ABC and evolutionary algorithm


## Experimental data
Three experimental datasets (contained in `../data`) under different temperatures are used in this section:
- `ExpGrowth.tsv` - the maximal specific growth rate in aerobic (Caspeta L., et al. Mbio, 2015) and anaerobic (Zakhartsev M., et al. J. Therm. Biol., 2015) batch cultivations

- `Chemostat_exp_data.txt` - fluxes of carbon dioxide (CO2), ethanol and glucose in chemostat cultivations (Postmus J., J. Biol. Chem., 2008)  

- `model_enzyme_params.csv` - File containing the enzyme parameters as determined by (Li, G *et al.*, Nature Communications, 2021). This file provides the basis for the unpermuted priors.

# Computational scripts

The computational scripts can be divided into three groups depending on their mode of operation:

- SLURM array jobs intended to be run on multiple nodes (with correspoinding SLURM script which most likely needs some modifications to run): `gem_smcabc_at_three_conditions_run.py` (`gem_smcabc_at_three_conditions_run.sh`), `gem_smcevo_tournament.py` (`gem_smcevo_tournament.sh`) and `gem_smcevo_truncation.py` (`gem_smcevo_truncation.sh`). Each of these scripts require their own preparation scripts `gem_smcabc_at_three_conditions_run.py`, `gem_smcevo_tournament_prepare.py` and `gem_smcevo_truncation_prepare.sh` to be run beforehand.
- Jobs which are indended to be run under SLURM, but where the script allows to be run the usual way on a single node: `evo_truncation_pca.py`, `evo_tournament_pca.py`, `gem_fva_at_three_conditions.py`, `reduced_pca.py` and `sample_pca.py`. The SLURM scripts have the same names just that they end in `.sh`.
- Jobs which are so small that no SLURM script is written, run them as usual scripts: `benchmark_performance.py`, `compute_particle_distances.py` and `reduce_data_size.py`.

Also note that the script must be run in a topologically sorted order to satisfy data dependencies. An attempt to illustrate the workflow is shown in `master_script.sh` which should in theory replicate the results, but it is written more for inspiration rather than actually being run. The script `benchmark_performance.py` is somewhat different and is described on its own section.

## SMC-ABC and evolutionary approach updates model parameters

These are the main scripts for fitting the parameters

* `gem_smcabc_at_three_conditions_prepare.py` - Prepare priors for the SMC-ABC. Output: `../results/permuted_smcabc_res/simulation_skeleton.pkl`
* `gem_smcabc_at_three_conditions_run.py` - SMC-ABC update parameter space with all three observed datasets. Requires: `../results/permuted_smcabc_res/simulation_skeleton.pkl`. Output: `../results/permuted_smcabc_res/smcabc_gem_three_conditions_updated_{prior}_{simulation}_save_all_particles.pkl`
* `gem_smcabc_tournament_prepare.py` - Prepare evolutionary simulations with tournament replacement for the evolutionary algorithm. Output: `../results/evo_tournament/simulation_skeleton.pkl`
* * `gem_smcabc_truncation_prepare.py` - Prepare evolutionary simulations with truncation replacement for the evolutionary algorithm. Output: `../results/evo_truncation/simulation_skeleton.pkl`
* `gem_smcabc_tournament.py` - Fits parameters with the evolutionary algorithm using tournament replacement. The `locality` parameter refers to the numbers of nearest neighbors to consider in the replacement. Requires: `../results/evo_tournament/simulation_skeleton.pkl`. Output: `../results/evo_tournament/smcevo_gem_tournament_{locality}_{simulation}.pkl`
* `gem_smcabc_truncation.py` - Fits parameters with the evolutionary algorithm using truncation replacement. The `num_elites` parameters refers to the number of the elites in the replacement, while the rest are picked on random. Requires: `../results/evo_truncation/simulation_skeleton.pkl`. Output: `../results/evo_truncation/smcevo_gem_truncation_{num_elites}_{simulation}.pkl`


Above scripts need be run on high-performance clusters, and may take a few days.

## PCA and t-SNE ordinations
* `sample_pca.py` - Create PCA ordination for all Bayesian simulations. Requires: `../results/permuted_smcabc_res/simulation_skeleton.pkl` and `../results/permuted_smcabc_res/smcabc_gem_three_conditions_updated_{prior}_{simulation}_save_all_particles.pkl`. Output: `../results/permuted_smcabc_res/particle_df.pkl`, `../results/permuted_smcabc_res/combined_particle_df.pkl` and `../results/permuted_smcabc_res/pca_full_ordination.pkl`.
* `reduced_pca.py` - Created PCA ordination for Bayesian unpermuted prior simulation 1 and 2 and Bayesian permuted prior 1 simulation 1 and 2. Requires: `../results/permuted_smcabc_res/combined_particle_df.pkl`. Output: `../results/permuted_smcabc_res/pca_reduced_ordination.pkl`.
* `evo_pca.py` - Creates PCA ordination for Bayesian unpermuted prior simulation 1 and the evolutionary simulation. Requires: `../results/permuted_smcabc_res/simulation_skeleton.pkl`, `../results/permuted_smcabc_res/smcabc_gem_three_conditions_updated_{origin}_{status}_save_all_particles.pkl` and `../results/smcevo_gem_three_conditions_save_all_particles_refined.pkl`. Output: `../results/evo_pca_full_ordination.pkl`, `../results/evo_combined_particle_df.pkl` and `../results/evo_combined_particle_df.pkl`.


## FVA analysis
* `gem_fva_at_three_conditions.py` - Run FVA on posterior particles from SMC-ABC. Requires: `../results/permuted_smcabc_res/simulation_skeleton.pkl` and `../results/permuted_smcabc_res/smcabc_gem_three_conditions_updated_{origin}_{status}_save_all_particles.pkl`. Output: `../results/permuted_smcabc_res/fva_at_three_conditions.pkl`.
* `evo_fva.py` - Run FVA on posterior particles from the evolutionary algorithm. Requires: `../results/smcevo_gem_three_conditions_save_all_particles_refined.pkl`. Output: `../results/evo_fva.pkl`.

## Summarize and visualize results

* `compute_particle_distances.py` - Computer weighted RMSD results for all particles. Requires: `../results/permuted_smcabc_res/combined_particle_df.pkl` and `../results/evo_particle_df.pkl`. Output: `../results/full_particle_RMSD.pkl`
* `reduce_data_size.py` - Utility script for stipping the computational results down to the bare minimum for creating visualizations. Requires: `../results/permuted_smcabc_res/simulation_skeleton.pkl` `../results/permuted_smcabc_res/smcabc_gem_three_conditions_updated_{origin}_{status}_save_all_particles.pkl`, `../results/evo_combined_particle_df.pkl`, `../results/permuted_smcabc_res/combined_particle_df.pkl`, `../results/smcevo_gem_three_conditions_save_all_particles_refined.pkl`, `../results/permuted_smcabc_res/fva_at_three_conditions.pkl` and `../results/evo_fva.pkl`.  Output: `../results/permuted_smcabc_res/combined_df_metadata.pkl`, `../results/evo_combined_df_metadata.pkl`, `../results/permuted_smcabc_res/distance_frame.pkl`, `../results/evo_distances.pkl` and `../results/aggregated_fva_res.pkl`
* `visualization_for_manuscript.ipynb` - A Jupyter notebook used for visualizing the results. Due to the bulk of the work offloaded to the other scripts this notebook should be light enough to use on any laptop or desktop computer. Requires: `../results/permuted_smcabc_res/distance_frame.pkl`, `../results/permuted_smcabc_res/combined_df_metadata.pkl`, `../results/permuted_smcabc_res/pca_full_ordination.pkl`, `../results/evo_combined_df_metadata.pkl`, `../results/evo_pca_full_ordination.pkl`, `../results/permuted_smcabc_res/pca_reduced_ordination.pkl`, `../results/evo_distances.pkl`, `../results/evo_tsne_res/tsne_skeleton.pkl`, `../results/evo_tsne_res/tsne_{i}.pkl`, `../results/full_particle_RMSD.pkl`, `../results/aggregated_fva_res.pkl`. Output: All plots used in the publication.

## Reproduce particle evaluation benchmarking

The script `benchmark_performance.py` benchmarks the performance of particle evaluation with reframed. Profiling with COBRApy is more tricky as it was used with an old version of the package and is now removed. However, fear not, the original particle evaluation framework is still in the git history in the branch `cobrapy_profile`. Run the script `benchmark_performance.sh` for automatically running the benchmarking procedure for both frameworks, it automatically handling swiching back and forth in git history.


