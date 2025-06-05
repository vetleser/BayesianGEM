#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import logging
from evo_etc import CrowdingDE
from sa_etc import SimulatedAnnealing

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

def build_a_dataframe_for_posterior_particles(file, r2_threshold = 0.98):
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



def perform_pca_on_parameters(df):
    epsilon = 1e-6

    logging.info(f"Shape of df is {df.shape}")

    # 1. normalize all columns to a standard normal distribution
    # df contains the trailing columns: r2, period, particle_ID, frame_ID, in total 4 columns to remove
    X = df.values[:,:-3]
    X_n = np.zeros_like(X)

    logging.info(f"Shape of X is {X.shape}")
    logging.info(f"Shape of X_n is {X_n.shape}")
    for i in range(X_n.shape[1]): X_n[:,i] = (X[:,i]-np.mean(X[:,i]))/(np.std(X[:,i]) + epsilon)
    pca = PCA(n_components=2)
    PCS = pca.fit_transform(X_n)
    logging.info(pca.explained_variance_ratio_)
    return PCS, pca.explained_variance_ratio_



outdir = '../results/sa'
#model_frame = load_pickle(f"{outdir}/sa_model_frame_R0985.pkl")

logging.info("Loading data")
#particle_dfs = list(map(build_a_dataframe_for_posterior_particles,model_frame.outfile))
logging.info("Augmenting data labeling")
df_1 = load_pickle(f"{outdir}/smcsa_gem_june2_0.1_0.5_0_df.pkl")
df_2 = load_pickle(f"{outdir}/smcsa_gem_june2_0.1_0.5_1_df.pkl")

logging.info(f"Dataframe 1 head: {df_1.head()}")
logging.info(f"Dataframe 2 head: {df_2.head()}")

treshold = 0.98
df_1 = df_1[df_1['r2'] > treshold]
df_2 = df_2[df_2['r2'] > treshold]
df_dict = {}
for i, particle_df in enumerate([df_1, df_2]):
    df_dict[i] = particle_df

logging.info("Combining dataframes")
combined_df = combine_dataframes(df_dict)
dump_pickle(combined_df,f"{outdir}/sa_combined_df_R098.pkl")
logging.info(f"Combined DataFrame shape: {combined_df.shape}")
logging.info("DataFrame head:")
logging.info(combined_df.head())    

logging.info("Performing PCA")
pca_ordination = perform_pca_on_parameters(combined_df)
dump_pickle(pca_ordination,f"{outdir}/sa_pca_R098.pkl")



logging.info("DONE")
