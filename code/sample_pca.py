#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import evo_etc
import multiprocessing

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

def build_a_dataframe_for_all_particles(results, n_priors = 128, r2_threshold = 0.9, convert_results=False):
    columns = list(results.all_particles[0].keys())
    columns.sort()
    print("Iterating over particles")
    data = list()
    for p in results.all_particles:
        if convert_results:
            p_refined = evo_etc.convert_to_raw_particle(p)
        else:
            p_refined = p
        data.append([p_refined[k] for k in columns])
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
        df["model"] = label
        df.reset_index()
    print("Labelling done")
    return pd.concat(augmented_df_dict.values(), ignore_index=True)

def perform_pca_on_parameters(df):
    # 1. normalize all columns to a standard normal distribution
    X, y, model, period = df.values[:,:-3], df.values[:,-3], df.values[:,-2], df.values[:,-1]
    X_n = np.zeros_like(X)    
    for i in range(X_n.shape[1]): X_n[:,i] = (X[:,i]-np.mean(X[:,i]))/np.std(X[:,i])
    pca = PCA(n_components=2)
    PCS = pca.fit_transform(X_n)
    print(pca.explained_variance_ratio_)
    return PCS, pca.explained_variance_ratio_

model_frame = load_pickle("../results/permuted_smcabc_res/result_model_frame.pkl")

def df_dict_helper(series):
        df_dict = {index: value for index, value in series.iteritems()}
        combined_df = combine_dataframes_for_models(df_dict)
        PCS, EV = perform_pca_on_parameters(combined_df)
        # This is a dirty trick to prevent variable unpacking and throwing away
        # the EV part of the result
        return PCS, EV


with multiprocessing.Pool(4) as p:
    job_dict = {name: group[name] for name, group in model_frame["particle_df"].groupby(level="origin")}
    raw_results = p.map(df_dict_helper, job_dict.values())
    pca_replication_ordinations = dict(zip(job_dict.keys(), raw_results))

dump_pickle(pca_replication_ordinations,"../results/pca_replication_ordinations.pkl")

original_particle_df_dict = {name: group[name]["original"] for name, group in model_frame["particle_df"].groupby(level="origin")}
combined_original_df = combine_dataframes_for_models(original_particle_df_dict)
pca_original_ordination = perform_pca_on_parameters(combined_original_df)
dump_pickle(pca_original_ordination,"../results/pca_original_ordination.pkl")
