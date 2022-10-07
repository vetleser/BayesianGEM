#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
import logging
import GEMS
from evo_etc import GA

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

def extract_distances_from_simulation(filename):
    sim_res = load_pickle(filename=filename)
    return sim_res.all_distances

def extract_distances_and_population_from_simulation(filename):
    Yobs_batch = GEMS.aerobic_exp_data()
    Yobs_chemo = GEMS.chemostat_exp_data()
    dfae_batch,dfan_batch =GEMS.load_exp_batch_data('../data/ExpGrowth.tsv')
    sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]
    Yobs_batch_an = {'data':dfan_batch.loc[sel_temp,'r_an'].values}
    Yobs = {'rae':Yobs_batch['data'],
            'ran':Yobs_batch_an['data']}
    sim_res: GA  = load_pickle(filename=filename)
    reduced_distances = [GEMS.distance_2(Yobs,res) for res in sim_res.all_simulated_data]
    return sim_res.all_distances,reduced_distances, sim_res.population

bayesian_simulation_skeleton = load_pickle("../results/permuted_smcabc_res/simulation_skeleton.pkl")

bayesian_simulation_skeleton["all_distances"] = list(map(extract_distances_from_simulation,bayesian_simulation_skeleton["outfile"]))
dump_pickle(bayesian_simulation_skeleton, "../results/permuted_smcabc_res/distance_frame.pkl")

evo_truncation_simulation_skeleton = load_pickle("../results/evo_truncation/simulation_skeleton.pkl")

evo_truncation_simulation_skeleton["all_distances"], evo_truncation_simulation_skeleton["reduced_distances"], evo_truncation_simulation_skeleton["population"] = zip(*list(map(extract_distances_and_population_from_simulation, evo_truncation_simulation_skeleton["outfile"])))

dump_pickle(evo_truncation_simulation_skeleton, "../results/evo_truncation/distance_frame.pkl")

evo_tournament_simulation_skeleton = load_pickle("../results/evo_tournament/simulation_skeleton.pkl")

evo_tournament_simulation_skeleton["all_distances"], evo_tournament_simulation_skeleton["reduced_distances"], evo_tournament_simulation_skeleton["population"] = zip(*list(map(extract_distances_and_population_from_simulation, evo_tournament_simulation_skeleton["outfile"])))

dump_pickle(evo_tournament_simulation_skeleton, "../results/evo_tournament/distance_frame.pkl")


bayesian_fva_results = load_pickle("../results/permuted_smcabc_res/fva_at_three_conditions.pkl")
flattened_df_list = []
signature_reactions_ids = list({'PDH': 'r_0961No1', 'FBA': 'r_0450No1', 'FCO': 'r_0438No1', 'PSP': 'r_0917No1', 'SHK': 'r_0997No1', 'GRW': 'r_2111'}.values())
for _, row in bayesian_fva_results.drop(columns=["particle"]).iterrows():
    raw_df = row["fva_res"]
    df = raw_df[np.isin(raw_df["reaction"],signature_reactions_ids)]
    df["origin"] = row["origin"]
    df["status"] = row["status"]
    flattened_df_list.append(df)

combined_fva_frame = (
    pd.concat(flattened_df_list).
    assign(range= lambda df: df["maximum"] - df["minimum"],
                                                         midpoint= lambda df: (df["maximum"] + df["minimum"]) / 2).
    drop(columns=["minimum", "maximum"])
        )

aggregated_fva_res = (
    combined_fva_frame.replace([np.inf, -np.inf],np.nan).
    dropna(how="all").
    groupby(["origin","status","condition","reaction","T"]).
    agg(["mean","min","max","std","count"])
                     )

dump_pickle(aggregated_fva_res,"../results/aggregated_fva_res.pkl")
