#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
import evo_etc
from sklearn.decomposition import PCA
import multiprocessing

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

def build_a_dataframe_for_all_particles(results, n_priors = 128, r2_threshold = 0.9):
    columns = list(results.all_particles[0].keys())
    columns.sort()
    print("Iterating over particles")
    data = list()
    for p in results.all_particles:
        data.append([p[k] for k in columns])
    print("Creating Data Frame")
    df = pd.DataFrame(data=data,columns=columns)
    df['r2'] = results.all_distances
    print(df.shape)
    
    # Remove samples with a R2 score smaller than -3
    print("Doing filtering and labelling of Data Frame")
    df['r2'] = -df['r2']
    sel_index = df.index[df['r2']>-3]    
    df = df.loc[sel_index,:]
    df["period"] = "Intermediate"
    df.loc[:n_priors,"period"] = "Prior"
    df.loc[df["r2"] > r2_threshold,"period"] = 'Posterior'
    print(df.shape)

    return df

def combine_dataframes_for_models(df_dict):
    # augmented_df_list =[ df.assign(model = lambda df: label)  for df, label in zip(df_list, index)]
    augmented_df_dict = {label: df.copy() for label, df in df_dict.items()}
    print("Copying done")
    for label, df in augmented_df_dict.items():
        df["origin"] = label[0]
        df["status"] = label[1]
        df.reset_index()
    print("Labelling done")
    return pd.concat(augmented_df_dict.values(), ignore_index=True)

def perform_pca_on_parameters(df):
    # 1. normalize all columns to a standard normal distribution
    X = df.values[:,:-3]
    X_n = np.zeros_like(X)    
    for i in range(X_n.shape[1]): X_n[:,i] = (X[:,i]-np.mean(X[:,i]))/np.std(X[:,i])
    pca = PCA(n_components=2)
    PCS = pca.fit_transform(X_n)
    print(pca.explained_variance_ratio_)
    return PCS, pca.explained_variance_ratio_

model_frame = load_pickle("../results/permuted_smcabc_res/result_model_frame.pkl")
with multiprocessing.Pool(10) as p:
    particle_df_map = p.map(build_a_dataframe_for_all_particles,model_frame.model)
    model_frame["particle_df"] = list(particle_df_map)

dump_pickle(model_frame["particle_df"], "../results/permuted_smcabc_res/particle_df.pkl")

full_particle_df_dict = {distribution: df for distribution, df in model_frame["particle_df"].iteritems()}
combined_df = combine_dataframes_for_models(full_particle_df_dict)
dump_pickle(combined_df, "../results/permuted_smcabc_res/combined_particle_df.pkl")
pca_ordination = perform_pca_on_parameters(combined_df)
dump_pickle(pca_ordination,"../results/permuted_smcabc_res/pca_full_ordination.pkl")
