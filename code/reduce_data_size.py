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
simulation_skeleton = load_pickle("../results/permuted_smcabc_res/simulation_skeleton.pkl")

def extract_distances_from_simulation(filename):
    sim_res = load_pickle(filename=filename)
    return sim_res.all_distances

simulation_skeleton["all_distances"] = list(map(extract_distances_from_simulation,simulation_skeleton["outfile"]))
dump_pickle(simulation_skeleton, "../results/permuted_smcabc_res/distance_frame.pkl")

evo_sim_res = load_pickle('../results/smcevo_gem_three_conditions_save_all_particles_refined.pkl')
dump_pickle(evo_sim_res.all_distances, "../results/evo_distances.pkl")