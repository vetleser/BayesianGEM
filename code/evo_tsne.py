#!/usr/bin/env python
# coding: utf-8

# This script is similar to evo_pca.py and has a direct dependence on it, but does t-SNE instead of PCA
# t-SNE is repeated on 12 different perplexities which results are intended to be compared afterwards. 
# Remember to prepare the array jobs by running evo_tsne_prepare.py first
import pickle
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import os
import logging


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

def perform_tsne_on_parameters(df, perplexity):
    # 1. normalize all columns to a standard normal distribution
    X = df.values[:,:-3]
    X_n = np.zeros_like(X)    
    for i in range(X_n.shape[1]): X_n[:,i] = (X[:,i]-np.mean(X[:,i]))/np.std(X[:,i])
    tsne = TSNE(2, perplexity=perplexity, learning_rate="auto", init="pca")
    return tsne.fit_transform(X_n)

task_idx = int(os.environ["SLURM_ARRAY_TASK_ID"])
logging.info("Loading data")
combined_df = load_pickle("../results/evo_combined_particle_df.pkl")
tsne_frame = load_pickle("../results/evo_tsne_res/tsne_skeleton.pkl")
logging.info(f"Task index is {task_idx}")
perplexity = tsne_frame["perplexity"][task_idx]
logging.info(f"Perplexity index is {perplexity}")
outfile = tsne_frame["outfile"][task_idx]
logging.info(f"Outfile is {outfile}")
logging.info("Running t-SNE")
ordination = perform_tsne_on_parameters(df=combined_df,perplexity=perplexity)
logging.info("Saving results")
dump_pickle(ordination,outfile)
logging.info("DONE")
