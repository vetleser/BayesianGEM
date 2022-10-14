#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# The purpose of this script is to trim the computational results from the simulations and subsequent PCA ordinations
# to the bare minimum required in the Jupyter notebook analyzing and visualizing the results.
# This is to make it fast and practical to run the notebook with laptop-grad hardware

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

bayesian_combined_df = load_pickle("../results/permuted_smcabc_res/combined_particle_df.pkl")
# We only use the period (Prior, Intermediate or Posterior), prior (unpermuted and permuted 0-2) and model (Simulation 1 or 2)
bayesian_combined_df_metadata = bayesian_combined_df[["period","origin","status"]]
dump_pickle(bayesian_combined_df_metadata,"../results/permuted_smcabc_res/combined_df_metadata.pkl")

reduced_bayesian_combined_df = load_pickle("../results/reduced_smcabc_res/combined_particle_df.pkl")
reduced_bayesian_combined_df_metadata = bayesian_combined_df[["period","simulation"]]
dump_pickle(reduced_bayesian_combined_df_metadata,"../results/reduced_smcabc_res/combined_df_metadata.pkl")

def extract_distances_from_simulation(filename):
    sim_res = load_pickle(filename=filename)
    return sim_res.all_distances

def extract_distances_and_population_from_simulation(filename):
    sim_res = load_pickle(filename=filename)
    return sim_res.all_distances, sim_res.population

bayesian_simulation_skeleton = load_pickle("../results/permuted_smcabc_res/simulation_skeleton.pkl")

bayesian_simulation_skeleton["all_distances"] = list(map(extract_distances_from_simulation,bayesian_simulation_skeleton["outfile"]))
dump_pickle(bayesian_simulation_skeleton, "../results/permuted_smcabc_res/distance_frame.pkl")

reduced_bayesian_simulation_skeleton = load_pickle("../results/reduced_smcabc_res/simulation_skeleton.pkl")
reduced_bayesian_simulation_skeleton["all_distances"] = list(map(extract_distances_from_simulation,reduced_bayesian_simulation_skeleton["outfile"]))
dump_pickle(reduced_bayesian_simulation_skeleton, "../results/reduced_smcabc_res/distance_frame.pkl")

evo_combined_df = load_pickle("../results/permuted_smcabc_res/combined_particle_df.pkl")
# We only use the period (Prior, Intermediate or Posterior), prior (unpermuted and permuted 0-2) and model (Simulation 1 or 2)
evo_combined_df_metadata = evo_combined_df[["particle_ID","period","method","frame_ID"]]
dump_pickle(evo_combined_df_metadata,"../results/evo_combined_df_metadata.pkl")

evo_truncation_simulation_skeleton = load_pickle("../results/evo_truncation/simulation_skeleton.pkl")

evo_truncation_simulation_skeleton["all_distances"], evo_truncation_simulation_skeleton["population"] = zip(*list(map(extract_distances_and_population_from_simulation, evo_truncation_simulation_skeleton["outfile"])))

dump_pickle(evo_truncation_simulation_skeleton, "../results/evo_truncation/distance_frame.pkl")

evo_tournament_simulation_skeleton = load_pickle("../results/evo_tournament/simulation_skeleton.pkl")

evo_tournament_simulation_skeleton["all_distances"], evo_tournament_simulation_skeleton["population"] = zip(*list(map(extract_distances_and_population_from_simulation, evo_tournament_simulation_skeleton["outfile"])))

dump_pickle(evo_tournament_simulation_skeleton, "../results/evo_tournament/distance_frame.pkl")

def aggregate_fva_results(result_df,simulation_attributes):
    flattened_df_list = []
    for _, row in result_df.drop(columns=["particle"]).iterrows():
        raw_df = row["fva_res"]
        df = raw_df
        for attribute in simulation_attributes:
            df[attribute] = row[attribute]
        flattened_df_list.append(df)

    combined_fva_frame = (
        pd.concat(flattened_df_list).
        assign(range= lambda df: df["maximum"] - df["minimum"],
                                                            midpoint= lambda df: (df["maximum"] + df["minimum"]) / 2).
        drop(columns=["minimum", "maximum"])
            )
    group_by_variables = simulation_attributes.extend(["condition","reaction","T"])
    aggregated_fva_res = (
        combined_fva_frame.replace([np.inf, -np.inf],np.nan).
        dropna(how="all").
        groupby(group_by_variables).
        agg(["mean","min","max","std","count"])
                        )
    return aggregated_fva_res

bayesian_fva_results = load_pickle("../results/permuted_smcabc_res/fva_at_three_conditions.pkl")
bayesian_aggregated_fva_results = aggregate_fva_results(bayesian_fva_results,["origin","status"])

reduced_bayesian_fva_results = load_pickle("../results/permuted_smcabc_res/fva_at_two_conditions.pkl")
reduced_bayesian_aggregated_fva_results = aggregate_fva_results(reduced_bayesian_fva_results,["simulation"])


dump_pickle(bayesian_aggregated_fva_results,"../results/aggregated_fva_res.pkl")

dump_pickle(reduced_bayesian_aggregated_fva_results,"../results/reduced_aggregated_fva_res.pkl")

evo_fva_results = load_pickle("../results/evo_fva.pkl")
truncation_fva_results = evo_fva_results["truncation"]
tournament_fva_results = evo_fva_results["tournament"]
truncation_aggregated_fva_results = aggregate_fva_results(truncation_fva_results,["num_elites","simulation"])
tournament_aggregated_fva_results = aggregate_fva_results(tournament_fva_results,["localty","simulation"])
dump_pickle(truncation_aggregated_fva_results, "../results/truncation_aggregated_fva_res.pkl")
dump_pickle(tournament_fva_results, "../results/tournament_aggregated_fva_res.pkl")
