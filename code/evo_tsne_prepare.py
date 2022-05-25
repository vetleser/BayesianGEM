#!/usr/bin/env python
# coding: utf-8

# This script is similar to evo_pca.py and has a direct dependence on it, but does t-SNE instead of PCA
# t-SNE is repeated on 12 different perplexities which results are intended to be compared afterwards. 
import pickle
import pandas as pd
import logging


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

perplexities = [5, 10, 30, 60, 100, 200, 300, 500, 1000, 2000, 5000, 10000]
outdir = '../results/evo_tsne_res'
tsne_frame = pd.DataFrame({"perplexity": perplexities})
tsne_frame["outfile"] = [f"{outdir}/tsne_{i}.pkl" for i, _ in enumerate(perplexities)]
dump_pickle(tsne_frame, f"{outdir}/tsne_skeleton.pkl")
logging.info("Setup DONE")
