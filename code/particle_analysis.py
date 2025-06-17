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

combined_normalized_importance = load_pickle(f"{outdir}/flux_analysis/combined_normalized_importance.pkl")



logging.info("Loading combined df")
#df = load_pickle(f"{outdir}/df_simulation_0_R098.pkl")
df = load_pickle(f"{outdir}/evo_combined_df_R098.pkl")


logging.info(f"Parameter data: \n{df}")

# Collect all columns that end with _Tm, _Topt, or _dCpt
tm_cols = [col for col in df.columns if col.endswith('_Tm')]
topt_cols = [col for col in df.columns if col.endswith('_Topt')]
dcpt_cols = [col for col in df.columns if col.endswith('_dCpt')]

# Compute global min and max for each parameter type
tm_min, tm_max = df[tm_cols].min().min(), df[tm_cols].max().max()
topt_min, topt_max = df[topt_cols].min().min(), df[topt_cols].max().max()
dcpt_min, dcpt_max = df[dcpt_cols].min().min(), df[dcpt_cols].max().max()

# Log or print the result
logging.info(f"Tm:    min = {tm_min:.2f}, max = {tm_max:.2f}")
logging.info(f"Topt:  min = {topt_min:.2f}, max = {topt_max:.2f}")
logging.info(f"dCpt:  min = {dcpt_min:.2f}, max = {dcpt_max:.2f}")

enzyme_id = 'P08566'
# Filter the DataFrame for the specific enzyme ID
df_enzyme = df[[col for col in df.columns if col.startswith(enzyme_id)]]
logging.info(f"Filtered DataFrame for enzyme {enzyme_id}: \n{df_enzyme}")
# Compute min and max for the specific enzyme
tm_cols = [col for col in df_enzyme.columns if col.endswith('_Tm')]
topt_cols = [col for col in df_enzyme.columns if col.endswith('_Topt')]
dcpt_cols = [col for col in df_enzyme.columns if col.endswith('_dCpt')]
tm_min_enzyme, tm_max_enzyme = df_enzyme[enzyme_id+"_Tm"].min().min(), df_enzyme[tm_cols].max().max()
topt_min_enzyme, topt_max_enzyme = df_enzyme[topt_cols].min().min(), df_enzyme[topt_cols].max().max()
dcpt_min_enzyme, dcpt_max_enzyme = df_enzyme[dcpt_cols].min().min(), df_enzyme[dcpt_cols].max().max()
# Log or print the result for the specific enzyme
logging.info(f"{enzyme_id} Tm:    min = {tm_min_enzyme:.2f}, max = {tm_max_enzyme:.2f}")
logging.info(f"{enzyme_id} Topt:  min = {topt_min_enzyme:.2f}, max = {topt_max_enzyme:.2f}")
logging.info(f"{enzyme_id} dCpt:  min = {dcpt_min_enzyme:.2f}, max = {dcpt_max_enzyme:.2f}")

df_sa = load_pickle(f"../results/sa/smcsa_gem_june2_0.1_0.5_0_df.pkl")
logging.info(f"SA df: \n{df_sa}")

tm_min_sa, tm_max_sa = df_sa[tm_cols].min().min(), df_sa[tm_cols].max().max()
topt_min_sa, topt_max_sa = df_sa[topt_cols].min().min(), df_sa[topt_cols].max().max()
dcpt_min_sa, dcpt_max_sa = df_sa[dcpt_cols].min().min(), df_sa[dcpt_cols].max().max()
# Log or print the result for SA
logging.info(f"SA Tm:    min = {tm_min_sa:.2f}, max = {tm_max_sa:.2f}")
logging.info(f"SA Topt:  min = {topt_min_sa:.2f}, max = {topt_max_sa:.2f}")
logging.info(f"SA dCpt:  min = {dcpt_min_sa:.2f}, max = {dcpt_max_sa:.2f}")

df_sa_098 = df_sa[df_sa['r2'] > 0.98]
logging.info(f"SA df with r2 > 0.98: \n{df_sa_098}")

tm_min_sa_098, tm_max_sa_098 = df_sa_098[tm_cols].min().min(), df_sa_098[tm_cols].max().max()
topt_min_sa_098, topt_max_sa_098 = df_sa_098[topt_cols].min().min(), df_sa_098[topt_cols].max().max()
dcpt_min_sa_098, dcpt_max_sa_098 = df_sa_098[dcpt_cols].min().min(), df_sa_098[dcpt_cols].max().max()
# Log or print the result for SA with r2 > 0.90
logging.info(f"SA Tm (r2 > 0.98):    min = {tm_min_sa_098:.2f}, max = {tm_max_sa_098:.2f}")
logging.info(f"SA Topt (r2 > 0.98):  min = {topt_min_sa_098:.2f}, max = {topt_max_sa_098:.2f}")
logging.info(f"SA dCpt (r2 > 0.98):  min = {dcpt_min_sa_098:.2f}, max = {dcpt_max_sa_098:.2f}")


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
    max_entropy = math.log2(len(bin_counts)) if len(bin_counts) > 1 else 1
    normalized_entropy = param_entropy / max_entropy
    logging.info(f"Entropy for {col} is normalized: {normalized_entropy}")
    logging.info(f"Length of bin_counts for {col} is: {len(bin_counts)}")

    # if param_entropy < 3.5:
    #     logging.info(f"Entropy for {col} is: {param_entropy}")
    x = [interval.mid for interval in bin_counts.index]
    y = bin_counts.values
    if param_entropy > 0:
        logging.info(f" Y.max() is: {y.max()} for {col}")
        logging.info(f"Entropy for {col} is: {param_entropy}")
        logging.info(f"plot_counter is: {plot_counter}")
        plot_counter += 1

        plt.figure(figsize=(12, 5))
        plt.bar(x, y, width=step, align='center')

        if col.endswith("_Tm"):
            plt.title(rf'Distribution of $T_m$ for Enzyme {col.split("_")[0]}', fontsize=16, fontweight='bold')
            plt.xlabel(r'$T_m$ (°K)')
        elif col.endswith("_Topt"):
            plt.title(rf'Distribution of $T_{{opt}}$ for Enzyme {col.split("_")[0]}', fontsize=16, fontweight='bold')
            plt.xlabel(r'$T_{{opt}}$ (°K)')
        else:
            plt.title(rf'Distribution of $dC_p\ddag$ for Enzyme {col.split("_")[0]}', fontsize=16, fontweight='bold')
            plt.xlabel(r'$\Delta C_p^\ddag$  (J/mol/K)')
        plt.ylabel('Count')
        #plt.title(f'Distribution of {col}. Entropy: {param_entropy:.4f}')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        plt.savefig(f"../figures/analysis/aaa_histogram_{col}_entropy.png")


# tot_entropy = {k: 0.0 for k in combined_normalized_importance.keys()}
# for key in df.columns.tolist():
#     if key.split('_')[0] not in protein_IDs:
#         logging.info(f"Skipping {key} as it is not in protein_IDs")
#         continue
#     start = math.floor(df[key].min())
#     end = math.ceil(df[key].max())
#     step = 0.5
#     if key.endswith("_dCpt"):
#         #logging.info("This is a dCpt column")
#         step = 100

#     # Create bin edges using numpy
#     bins = np.arange(start, end + step, step)

#     # Use pd.cut with these bins
#     binned = pd.cut(df[key], bins=bins)
#     bin_counts = binned.value_counts().sort_index()
#     param_entropy = entropy(bin_counts)
#     max_entropy = math.log2(len(bin_counts)) if len(bin_counts) > 1 else 1
#     normalized_entropy = param_entropy / max_entropy
#     tot_entropy[key.split('_')[0]] += normalized_entropy

# bar_colors = []
# for k in tot_entropy.keys():
#     importance = combined_normalized_importance.get(k, None)
#     if importance == 0:
#         bar_colors.append('red')     # Color for importance = 0
#     elif importance == 1:
#         bar_colors.append('green')   # Color for importance = 1
#     else:
#         bar_colors.append('blue')    # Default color for in-between values

# entropy_0 = []
# entropy_1 = []
# entropy_mid = []

# for k, entropy_val in tot_entropy.items():
#     importance = combined_normalized_importance.get(k, None)
#     if importance is None:
#         continue
#     if importance == 0:
#         entropy_0.append(entropy_val)
#     elif importance == 1:
#         entropy_1.append(entropy_val)
#     elif 0 < importance < 1:
#         entropy_mid.append(entropy_val)

# # Compute averages
# avg_entropy_0 = np.mean(entropy_0) if entropy_0 else float('nan')
# avg_entropy_1 = np.mean(entropy_1) if entropy_1 else float('nan')
# avg_entropy_mid = np.mean(entropy_mid) if entropy_mid else float('nan')
# std_entropy_0 = np.std(entropy_0) if entropy_0 else float('nan')
# std_entropy_1 = np.std(entropy_1) if entropy_1 else float('nan')
# std_entropy_mid = np.std(entropy_mid) if entropy_mid else float('nan')
# logging.info(f"Average entropy for importance = 0: {avg_entropy_0} (std: {std_entropy_0})")
# logging.info(f"Average entropy for importance = 1: {avg_entropy_1} (std: {std_entropy_1})")
# logging.info(f"Average entropy for importance between 0 and 1: {avg_entropy_mid} (std: {std_entropy_mid})")



# logging.info(f"Length of tot_entropy: {len(tot_entropy)}")

#logging.info(f"Total entropy for each parameter: {tot_entropy}")
plot_counter = 0

# cols1 = ['P08566', 'Q99190', 'P38286', 'P40857', 'P47176', 'P00815', 'P05375', 'P07245', 'P40319', 'P36010']

cols1 = [enzyme_id]

for prot in cols1:
    param = prot + '_Tm'
    plot_histogram_entropy(param)
    param = prot + '_Topt'
    plot_histogram_entropy(param)
    param = prot + '_dCpt'
    plot_histogram_entropy(param)
    plt.close()

# cols = df.columns.tolist()
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