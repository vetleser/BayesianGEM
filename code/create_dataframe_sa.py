#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import logging
from evo_etc import CrowdingDE
from sa_etc import SimulatedAnnealing
import os

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

def build_a_dataframe_for_posterior_particles(file, r2_threshold = 0.90):
    logging.info(f"Loading file: {file}")
    results: SimulatedAnnealing = load_pickle(file)
    columns = list(results.all_particles[0].keys())
    columns.sort()
    logging.info("Iterating over particles")
    data = list()
    for p in results.all_particles:
        data.append([p[k] for k in columns])
    logging.info("Creating Data Frame")
    df = pd.DataFrame(data=data,columns=columns)
    df['r2'] = results.all_distances
    # Running number assigned to particles to keep track of them when comparing with original data
    df['particle_ID'] = list(range(len(results.all_particles)))
    logging.info(df.shape)
    
    logging.info("Doing filtering and labelling of Data Frame")
    # We need to negate the results due to the fact that they
    # are orignally taken to mean distances which are to be minimized
    df['r2'] = -df['r2']
    df = df[df['r2'] > r2_threshold]
    logging.info(df.shape)
    return df

def combine_dataframes(df_dict):
    for label, df in df_dict.items():
        df["frame_ID"] = label
    return pd.concat(df_dict.values(),ignore_index=True)


outdir = '../results/sa'
file = f'{outdir}/sa_simulation_skeleton.pkl'
logging.info(f"File size: {os.path.getsize(file)} bytes")

model_frame: pd.DataFrame = load_pickle(file)
model_frame["frame_ID"] = range(model_frame.shape[0])
# Load all from the model frame
reduced_model_frame = (
    model_frame.
    set_index(["end_exploration"]).
    loc[(0.25, 0.5, 0.75), :]
)
logging.info(f"Reduced model frame shape: {reduced_model_frame.shape}")
logging.info(f"Reduced model frame columns: {reduced_model_frame.columns}")
logging.info(f"Reduced model frame:\n{reduced_model_frame}")
# logging.info("Loading data")
# particle_dfs = list(map(build_a_dataframe_for_posterior_particles,reduced_model_frame.outfile))
# logging.info("Augmenting data labeling")
# df_dict = {}
# for i, particle_df in zip(reduced_model_frame["frame_ID"],particle_dfs):
#     df_dict[i] = particle_df

# logging.info("Combining dataframes")
# combined_df = combine_dataframes(df_dict)
# logging.info(f"Combined DataFrame shape: {combined_df.shape}")
# logging.info(f"Combined DataFrame columns: {combined_df}")
# dump_pickle(combined_df,f"{outdir}/sa_combined_df_R090_final.pkl")

combined_df = load_pickle(f"{outdir}/sa_combined_df_R090_final.pkl")
logging.info(f"Combined DataFrame shape: {combined_df.shape}")
logging.info(f"Combined DataFrame columns: {combined_df.columns}")
r2_threshold = 0.97
logging.info(f"Filtering DataFrame with r2 threshold: {r2_threshold}")
filtered_df = combined_df[combined_df['r2'] > r2_threshold]
logging.info(f"Filtered DataFrame shape: {filtered_df.shape}")
logging.info(f"Filtered DataFrame columns: {filtered_df.columns}")
for i in range(0, 6):
    logging.info(f"Filtered DataFrame for frame_ID {i}:")
    logging.info(filtered_df[filtered_df['frame_ID'] == i].shape)
