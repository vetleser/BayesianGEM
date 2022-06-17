#!/usr/bin/env python
# coding: utf-8

import pickle
import logging
import pandas as pd
import numpy as np

epsilon = 1e-5
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))


bayesian_particle_df = load_pickle("../results/permuted_smcabc_res/combined_particle_df.pkl")
evolutionary_particle_df = load_pickle("../results/evo_particle_df.pkl")
evolutionary_particle_df["origin"] = "unpermuted"
evolutionary_particle_df["status"] = "evolutionary"
all_particle_df = (
    pd.concat((bayesian_particle_df,evolutionary_particle_df),ignore_index=True).
    query("period == 'Posterior'").
    drop(columns=["r2","period"]).
    assign(particle=lambda df: range(df.shape[0]))
    )
long_particle_df = all_particle_df.melt(id_vars=["origin","status","particle"],var_name="reaction",value_name="value")
all_standard_deviations = long_particle_df.drop(columns=["origin","status","particle"]).groupby("reaction").agg(np.std)["value"]
# To guard agains division by zero problems
all_standard_deviations[all_standard_deviations < epsilon] = epsilon
result = (
    long_particle_df.set_index("particle").groupby(["origin","status","reaction"]).
    apply(lambda df: ((df["value"] - df["value"].mean()) / all_standard_deviations[df.name[2]])**2).
    reset_index(drop=False).
    groupby(["origin","status", "particle"]).
    apply(lambda df: np.sqrt(df["value"].mean())).
    reset_index(drop=False,name="RMSD").
    drop(columns="particle")
)
dump_pickle(result, "../results/full_particle_RMSD.pkl")
