#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import logging
from evo_etc import GA

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

def build_a_dataframe_for_all_particles(file, r2_threshold = 0.9):
    results: GA = load_pickle(file)
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
    df['ID'] = list(range(len(results.all_particles)))
    logging.info(df.shape)
    
    logging.info("Doing filtering and labelling of Data Frame")
    # We need to negate the results due to the fact that they 
    # are orignally taken to mean distances which are to be minimized
    df['r2'] = -df['r2']
    df["period"] = "Intermediate"
    df.loc[np.array(results.birth_generation) == 0,"period"] = "Prior"
    df.loc[df["r2"] > r2_threshold,"period"] = 'Posterior'
    # Remove samples with a R2 score smaller than -3
    sel_index = df.index[df['r2']>-3]    
    df = df.loc[sel_index,:]
    logging.info(df.shape)

    return df

def perform_pca_on_parameters(df):
    epsilon = 1e-6
    # 1. normalize all columns to a standard normal distribution
    X = df.values[:,:-3]
    X_n = np.zeros_like(X)    
    for i in range(X_n.shape[1]): X_n[:,i] = (X[:,i]-np.mean(X[:,i]))/(np.std(X[:,i]) + epsilon)
    pca = PCA(n_components=2)
    PCS = pca.fit_transform(X_n)
    logging.info(pca.explained_variance_ratio_)
    return PCS, pca.explained_variance_ratio_



outdir = '../results/evo_truncation'
model_frame = load_pickle(f"{outdir}/simulation_skeleton.pkl")

model_frame.set_index(["num_elites","simulation"], inplace=True)


particle_df_map = list(map(build_a_dataframe_for_all_particles,model_frame.outfile))
model_frame["particle_df"] = particle_df_map
dump_pickle(model_frame["particle_df"], f"{outdir}/particle_df.pkl")
logging.info("Performing PCA")
pca_ordination = list(map(perform_pca_on_parameters,particle_df_map))
model_frame["pca_ordination"] = pca_ordination
dump_pickle(pca_ordination,f"{outdir}/pca_all_ordination.pkl")



logging.info("DONE")
