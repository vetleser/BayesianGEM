## Description
This folder contains the scripts that carry out all analysis reported in the paper. Some scripts assumes the job is run under the SLURM workload manager. While these scripts will require some modification to run under your perferred computing environment, such modifications should be relatively feasible to carry out.

## Backend scripts
* `etcpy` - the scripts for incorporating temperature parameters into the enzyme constrained Yeast 7.6 model, the detail introduction can be found from `etcpy/README.md`.
* `GEMS.py` - the functionality for simulating the aerobic, anaerobic growth rate in batch cultivation and aerobic fluxes in chemostat cultivation, as well as a list of distance functions used for SMC-ABC and the evolutionary approach. It also includes corresponding FVA functionality for the batch and chemostat cultivations.
* `abc_etc.py` - the functionality to perform SMC-ABC approach.
* `evo_etc` - Like `abc_etc.py`, but uses the CrowdingDE algorithm instead of SMC-ABC
* `permute_parameters.py` - Utility script for shuffling parameters and thus creating permuted priors
* `random_sampler.py` - Utility script for creating random numbers for SMC-ABC and evolutionary algorithm


## Experimental data
Three experimental datasets (contained in `../data`) under different temperatures are used in this section:
- `ExpGrowth.tsv` - the maximal specific growth rate in aerobic (Caspeta L., et al. Mbio, 2015) and anaerobic (Zakhartsev M., et al. J. Therm. Biol., 2015) batch cultivations

- `Chemostat_exp_data.txt` - fluxes of carbon dioxide (CO2), ethanol and glucose in chemostat cultivations (Postmus J., J. Biol. Chem., 2008)  

- `model_enzyme_params.csv` - File containing the enzyme parameters as determined by (Li, G *et al.*, Nature Communications, 2021). This file provides the basis for the unpermuted priors.

# Computational scripts

The computational scripts can be divided into three groups depending on their mode of operation:

- SLURM array jobs intended to be run on multiple nodes (with correspoinding SLURM script which most likely needs some modifications to run): `gem_smcabc_at_three_conditions_run.py` (`gem_smcabc_at_three_conditions_run.sh`), `gem_smcabc_at_two_conditions.py` (`gem_smcabc_at_two_conditions.sh`), and `gem_smcevo_run.py` (`gem_smcevo_run.sh`). Each of these scripts require their own preparation scripts `gem_smcabc_at_three_conditions_prepare.py`, `gem_smcabc_at_two_conditions_prepare.py`, and `gem_smcevo_prepare.py` to be run beforehand.
- Jobs which are indended to be run under SLURM, but where the script allows to be run the usual way on a single node: `evo_pca.py`, `gem_fva_at_three_conditions.py`. `gem_fva_at_two_conditions.py`, `evo_pca.py` and `sample_pca.py`. The SLURM scripts have the same names just that they end in `.sh`.
- Jobs which are so small that no SLURM script is written, run them as usual scripts: `toy_example.py`, `benchmark_performance.py` and `reduce_data_size.py`.

Also note that the script must be run in a topologically sorted order to satisfy data dependencies. An attempt to illustrate the workflow is shown in `master_script.sh` which should in theory replicate the results, but it is written more for inspiration rather than actually being run. The script `benchmark_performance.py` is somewhat different and is described on its own section.

## SMC-ABC and evolutionary approach updates model parameters

These are the main scripts for fitting the parameters

* `gem_smcabc_at_three_conditions_prepare.py` - Prepare priors for the SMC-ABC. Output: `../results/permuted_smcabc_res/simulation_skeleton.pkl`
* `gem_smcabc_at_two_conditions_prepare.py` - Prepare priors for the SMC-ABC for two condition without permuted priors. Output: `../results/reduced_smcabc_res/simulation_skeleton.pkl`
* `gem_smcabc_at_three_conditions_run.py` - Run Bayesian calculation method with all three observed datasets. Requires: `../results/permuted_smcabc_res/simulation_skeleton.pkl`. Output: `../results/permuted_smcabc_res/smcabc_gem_three_conditions_updated_{prior}_{simulation}_save_all_particles.pkl`
* `gem_smcabc_at_two_conditions.py` - Run Bayesian calculation method with the aerobic and anaerobic datasets. Requires: `../results/reduced_smcabc_res/simulation_skeleton.pkl`. Output: `../results/permuted_smcabc_res/smcabc_gem_two_conditions_{simulation}_save_all_particles.pkl`
* `gem_smcevo_prepare.py` - Prepare evolutionary simulations with the CrowdingDE algorithm. Output: `../results/crowdingDE/simulation_skeleton.pkl`
* `gem_smcevo_run.py` - Fits parameters with the evolutionary algorithm with the CrowdingDE algorithm. The `scaling_factor` parameter refers to the effect of the secondary parents, whereas the `crossover_prob` parameter determines the number of crossover with the two secondary parents.. Requires: `../results/crowdingDE/simulation_skeleton.pkl`. Output: `../results/crowdingDE/smcevo_gem_{scaling_factor}_{crossover_prob}_{simulation}.pkl`


Above scripts need be run on high-performance clusters, and may take a few days.

## PCA ordinations
* `sample_pca.py` - Create PCA ordination for Bayesian calculation method for all three datasets. Requires: `../results/permuted_smcabc_res/simulation_skeleton.pkl` and `../results/permuted_smcabc_res/smcabc_gem_three_conditions_updated_{prior}_{simulation}_save_all_particles.pkl`. Output: `../results/permuted_smcabc_res/particle_df.pkl`, `../results/permuted_smcabc_res/combined_particle_df.pkl`, `../results/permuted_smcabc_res/pca_reduced_ordination.pkl`, and `../results/permuted_smcabc_res/pca_full_ordination.pkl`.
* `reduced_pca.py` - Created PCA ordination for Bayesian calculation method excluding the chemostat dataset. Requires: `../results/reduced_smcabc_res/simulation_skeleton.pkl` and `../results/permuted_smcabc_res/smcabc_gem_two_conditions_{simulation}_save_all_particles.pkl`. Output: `../results/reduced_smcabc_res/pca_full_ordination.pkl` and `../results/reduced_smcabc_res/combined_particle_df.pkl`.
* `evo_pca.py` - Creates PCA ordination for CrowdingDE. Requires: `../results/crowdingDE/smcevo_gem_{scaling_factor}_{crossover_prob}_{simulation}.pkl`, and `../results/crowdingDE/simulation_skeleton.pkl`. Output: `../results/crowdingDE/evo_combined_df.pkl` and `../results/crowdingDE/evo_pca.pkl`.


## FVA analysis
* `gem_fva_at_three_conditions.py` - Run FVA on posterior particles from Bayesian calculation method. Requires: `../results/permuted_smcabc_res/simulation_skeleton.pkl` and `../results/permuted_smcabc_res/smcabc_gem_three_conditions_updated_{origin}_{status}_save_all_particles.pkl`. Output: `../results/permuted_smcabc_res/fva_at_three_conditions.pkl`.
* `gem_fva_at_two_conditions.py` - Run FVA on posterior particles from Bayesian calculation method excluding the chemostat dataset. Requires: `../results/reduced_smcabc_res/simulation_skeleton.pkl` and `../results/permuted_smcabc_res/smcabc_gem_two_conditions_{simulation}_save_all_particles.pkl`. Output: `../results/reduced_smcabc_res/fva_at_two_conditions.pkl`.
* `evo_fva.py` - Run FVA on posterior particles from the evolutionary algorithm. Requires: `../results/crowdingDE/smcevo_gem_{scaling_factor}_{crossover_prob}_{simulation}.pkl`, and `../results/crowdingDE/simulation_skeleton.pkl`.. Output: `../results/crowdingDE/evo_fva.pkl`.

## Summarize and visualize results

* `toy_example.py` - Script for running the toy example which compares the Bayesian calculation method and the evolutionary algorithm.
* `reduce_data_size.py` - Utility script for stipping the computational results down to the bare minimum for creating visualizations. Requires: `../results/permuted_smcabc_res/simulation_skeleton.pkl` `../results/permuted_smcabc_res/smcabc_gem_three_conditions_updated_{origin}_{status}_save_all_particles.pkl`, `../results/permuted_smcabc_res/combined_particle_df.pkl`,`../results/reduced_smcabc_res/simulation_skeleton.pkl`, `../results/permuted_smcabc_res/smcabc_gem_two_conditions_{simulation}_save_all_particles.pkl`, `../results/reduced_smcabc_res/combined_particle_df.pkl`, `../results/crowdingDE/smcevo_gem_{scaling_factor}_{crossover_prob}_{simulation}.pkl`, `../results/crowdingDE/simulation_skeleton.pkl`, `../results/crowdingDE/evo_combined_df.pkl`, `../results/permuted_smcabc_res/fva_at_three_conditions.pkl`, `../results/reduced_smcabc_res/fva_at_two_conditions.pkl`, and `../results/crowdingDE/evo_fva.pkl`.  Output: `../results/permuted_smcabc_res/combined_df_metadata.pkl`, `../results/reduced_smcabc_res/combined_df_metadata.pkl`, `../results/crowdingDE/evo_combined_df_metadata.pkl`, `../results/permuted_smcabc_res/distance_frame.pkl`, `../results/reduced_smcabc_res/distance_frame.pkl`, `../results/crowdingDE/distance_frame.pkl`, `../results/crowdingDE/evo_aggregated_fva_res.pkl`, `../results/aggregated_fva_res.pkl` and `../results/reduced_aggregated_fva_res.pkl`.
* `visualization_for_manuscript.ipynb` - A Jupyter notebook used for visualizing the results. Due to the bulk of the work offloaded to the other scripts this notebook should be light enough to use on any laptop or desktop computer. Requires: `../results/permuted_smcabc_res/combined_df_metadata.pkl`, `../results/reduced_smcabc_res/combined_df_metadata.pkl`, `../results/pca.pkl`, `../results/permuted_smcabc_res/distance_frame.pkl`, `../results/reduced_smcabc_res/distance_frame.pkl`, `../results/crowdingDE/distance_frame.pkl`, `../results/crowdingDE/evo_aggregated_fva_res.pkl`, `../results/aggregated_fva_res.pkl`, `../results/reduced_aggregated_fva_res.pkl`, ``../results/permuted_smcabc_res/pca_full_ordination.pkl`, `../results/reduced_smcabc_res/pca_full_ordination.pkl`, `../results/permuted_smcabc_res/pca_reduced_ordination.pkl`, `../results/crowdingDE/evo_combined_df.pkl`, and `../results/crowdingDE/evo_pca.pkl`. Output: All plots used in the publication.

## Reproduce particle evaluation benchmarking

The script `benchmark_performance.py` benchmarks the performance of particle evaluation with reframed. Profiling with COBRApy is more tricky as it was used with an old version of the package and is now removed. However, fear not, the original particle evaluation framework is still in the git history in the branch `cobrapy_profile`. Run the script `benchmark_performance.sh` for automatically running the benchmarking procedure for both frameworks, it automatically handling swiching back and forth in git history.
