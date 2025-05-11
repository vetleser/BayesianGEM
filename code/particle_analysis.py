#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
import logging
import evo_etc as CrowdingDE
import math

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D



logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)

#Files to look at:
#evo_pca, reduce_data_size

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
logging.info("BEGIN")
outdir = "../results/analysis"


# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))


logging.info("Loading protein IDs")
protein_IDs = pd.read_csv("../data/model_enzyme_params.csv").iloc[:,0].tolist()



logging.info("Loading combined df")
df = load_pickle(f"{outdir}/df_simulation_0_R098.pkl")

logging.info(f"Parameter data: \n{df}")



def plot_histogram(col):
    global plot_counter
    # Define bins (e.g., from min to max with step size 5)
    start = math.floor(df[col].min())
    end = math.ceil(df[col].max())
    step = 0.5
    if col.endswith("_dCpt"):
        #logging.info("This is a dCpt column")
        step = 100

    # Create bin edges using numpy
    bins = np.arange(start, end + step, step)

    # Use pd.cut with these bins
    binned = pd.cut(df[col], bins=bins)
    bin_counts = binned.value_counts().sort_index()
    x = [interval.mid for interval in bin_counts.index]
    y = bin_counts.values
    if y.max() > 500:
        logging.info(f" Y.max() is: {y.max()} for {col}")
        logging.info(f"plot_counter is: {plot_counter}")
        plot_counter += 1

        plt.figure(figsize=(12, 5))
        plt.bar(x, y, width=step, align='center')

        plt.xlabel('Temperature (°K)')
        plt.ylabel('Count')
        plt.title(f'Distribution of {col}')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        plt.savefig(f"../figures/analysis/histogram_{col}.png")
    if plot_counter ==10:
        logging.info("Too many bins, skipping histogram")
        return

col = 'O13525_Tm'

plot_counter = 0


cols = df.columns.tolist()
for col in cols:
    plot_histogram(col)
    plt.close()
    if plot_counter == 10:
        logging.info("Too many bins, skipping histogram")
        break