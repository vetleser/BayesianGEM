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

full_particle_df_dict = {distribution: df for distribution, df in model_frame["particle_df"].iteritems()}
combined_df = combine_dataframes_for_models(full_particle_df_dict)
dump_pickle(combined_df, "../results/permuted_smcabc_res/combined_particle_df.pkl")
pca_ordination = perform_pca_on_parameters(combined_df)
dump_pickle(pca_ordination,"../results/permuted_smcabc_res/pca_full_ordination.pkl")
