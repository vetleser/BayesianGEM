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

combined_df = load_pickle("../results/permuted_smcabc_res/combined_particle_df.pkl")
# We only use the period (Prior, Intermediate or Posterior), prior (unpermuted and permuted 0-2) and model (Simulation 1 or 2)
combined_df_metadata = combined_df[["period","origin","status"]]
dump_pickle(combined_df_metadata,"../results/permuted_smcabc_res/combined_df_metadata.pkl")

evo_combined_df = load_pickle("../results/evo_combined_particle_df.pkl")
evo_combined_df_metadata = evo_combined_df[["model","period"]]
dump_pickle(evo_combined_df_metadata, "../results/evo_combined_df_metadata.pkl")


bayesian_simulation_skeleton = load_pickle("../results/permuted_smcabc_res/simulation_skeleton.pkl")

def extract_distances_from_simulation(filename):
    sim_res = load_pickle(filename=filename)
    return sim_res.all_distances

bayesian_simulation_skeleton["all_distances"] = list(map(extract_distances_from_simulation,bayesian_simulation_skeleton["outfile"]))
dump_pickle(bayesian_simulation_skeleton, "../results/permuted_smcabc_res/distance_frame.pkl")

evo_simulation_skeleton = load_pickle("../results/permuted_smcevo_res/simulation_skeleton.pkl")

evo_simulation_skeleton["all_distances"] = list(map(extract_distances_from_simulation,evo_simulation_skeleton["outfile"]))

dump_pickle(evo_simulation_skeleton, "../results/permuted_smcevo_res/distance_frame.pkl")
evo_sim_res = load_pickle('../results/smcevo_gem_three_conditions_save_all_particles_refined.pkl')

bayesian_fva_results = load_pickle("../results/permuted_smcabc_res/fva_at_three_conditions.pkl")
# evo_fva_results = load_pickle("../results/evo_fva.pkl")
flattened_df_list = []
for _, row in bayesian_fva_results.drop(columns=["particle"]).iterrows():
    df = row["fva_res"].copy()
    df[""] = row["origin"]
    df["status"] = row["status"]
    flattened_df_list.append(df)

aggregated_fva_res = (bayesian_fva_results.assign(range= lambda df: df["maximum"] - df["minimum"],
                                                         midpoint= lambda df: (df["maximum"] + df["minimum"]) / 2).
    drop(columns=["minimum", "maximum"]).
    replace([np.inf, -np.inf],np.nan).
    dropna(how="all").
    groupby(["origin","status","condition","reaction","T"]).
    agg(["mean","min","max","std","count"])
        )


dump_pickle(aggregated_fva_res,"../results/aggregated_fva_res.pkl")
