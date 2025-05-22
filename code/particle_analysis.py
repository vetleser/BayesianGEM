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
from scipy.stats import entropy, spearmanr



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
#df = load_pickle(f"{outdir}/df_simulation_0_R098.pkl")
df = load_pickle(f"{outdir}/evo_combined_df_R098.pkl")


logging.info(f"Parameter data: \n{df}")



def plot_histogram_ymax(col):
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
    param_entropy = entropy(bin_counts)
    #logging.info(f"Entropy for {col} is: {param_entropy}")
    x = [interval.mid for interval in bin_counts.index]
    y = bin_counts.values
    if y.max() > 1000:
        logging.info(f" Y.max() is: {y.max()} for {col}")
        logging.info(f"Entropy for {col} is: {param_entropy}")
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

def plot_histogram_entropy(col):
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
    param_entropy = entropy(bin_counts)
    if param_entropy < 3.5:
        logging.info(f"Entropy for {col} is: {param_entropy}")
    x = [interval.mid for interval in bin_counts.index]
    y = bin_counts.values
    if param_entropy < 10:
        logging.info(f" Y.max() is: {y.max()} for {col}")
        logging.info(f"Entropy for {col} is: {param_entropy}")
        logging.info(f"plot_counter is: {plot_counter}")
        plot_counter += 1

        plt.figure(figsize=(12, 5))
        plt.bar(x, y, width=step, align='center')

        plt.xlabel('Temperature (°K)')
        plt.ylabel('Count')
        plt.title(f'Distribution of {col}. Entropy: {param_entropy:.4f}')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        plt.savefig(f"../figures/analysis/aa_histogram_{col}_entropy.png")

col = 'O13525_Tm'

plot_counter = 0

cols1 = ['P08566', 'Q99190', 'P38286', 'P40857', 'P47176', 'P00815', 'P05375', 'P07245', 'P40319', 'P36010']

for prot in cols1:
    param = prot + '_Tm'
    plot_histogram_entropy(param)
    plt.close()

cols = df.columns.tolist()
# for col in cols:
#     plot_histogram_entropy(col)
#     plt.close()
#     if plot_counter == 10:
#         logging.info("Too many bins, skipping histogram")
#         break

# for col in cols:
#     if col.endswith("_Topt"):
#         plot_histogram_entropy(col)
#         plt.close()
        

#         if plot_counter == 10:
#             logging.info("Too many bins, skipping histogram")
#             break

# entropy_list = []
# for col in cols:
#     # if not (col.endswith("_Tm") or col.endswith("_Topt")):
#         # continue
#     # if not (col.endswith("_dCpt")):
#     #     continue
#     start = math.floor(df[col].min())
#     end = math.ceil(df[col].max())
#     step = 0.5
#     if col.endswith("_dCpt"):
#         #logging.info("This is a dCpt column")
#         step = 100

#     # Create bin edges using numpy
#     bins = np.arange(start, end + step, step)
#     logging.info(f"Number of bins for {col} are: {len(bins)}")

#     # Use pd.cut with these bins
#     binned = pd.cut(df[col], bins=bins)
#     bin_counts = binned.value_counts().sort_index()
#     param_entropy = entropy(bin_counts)
#     entropy_list.append(param_entropy)

# # Plotting the entropy values
# # Pair column names with their entropies
# #entropies = pd.Series(entropy_list, index=cols)

# logging.info(f"Entropy list minimum: {min(entropy_list)}")
# logging.info(f"Entropy list maximum: {max(entropy_list)}")

# plt.figure(figsize=(8,4))
# plt.bar(range(len(entropy_list)), entropy_list)
# plt.ylabel("Shannon entropy")
# plt.title("Entropy of each feature")
# #plt.xticks(rotation=45, ha="right")
# plt.tight_layout()
# plt.show()
# plt.savefig(f"../figures/analysis/entropy.png")

# logging.info("DONE")