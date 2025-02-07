#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import numpy as np
import logging
import evo_etc as CrowdingDE

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D



logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)

#Files to look at:
#evo_pca, reduce_data_size

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
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
file = load_pickle(f"{outdir}/evo_combined_df_R098.pkl")




logging.info("Creating DataFrames for std and mean")
std_df = file.iloc[:, :-3].std().to_frame().T  # Compute std and convert to a single-row DataFrame
mean_df = file.iloc[:, :-3].mean().to_frame().T  # Compute mean and convert to a single-row DataFrame

logging.info("Creating DataFrame with Protein-columns, each parameter STD as a row")
df_protein_std = pd.DataFrame(columns=protein_IDs)

for protein in protein_IDs:
    df_protein_std.loc[0, protein] = std_df.get(f"{protein}_Tm", pd.Series([None])).iloc[0]
    df_protein_std.loc[1, protein] = std_df.get(f"{protein}_Topt", pd.Series([None])).iloc[0]
    df_protein_std.loc[2, protein] = std_df.get(f"{protein}_dCpt", pd.Series([None])).iloc[0]

df_protein_std["Value_type"] = ["Tm", "T_opt", "dCpt"]

dump_pickle(df_protein_std, f"{outdir}/df_protein_std.pkl")
#print(df_protein_std)

logging.info("Creating DataFrame with Protein-columns, each parameter MEAN as a row")
df_protein_mean = pd.DataFrame(columns=protein_IDs)

for protein in protein_IDs:
    df_protein_mean.loc[0, protein] = mean_df.get(f"{protein}_Tm", pd.Series([None])).iloc[0]
    df_protein_mean.loc[1, protein] = mean_df.get(f"{protein}_Topt", pd.Series([None])).iloc[0]
    df_protein_mean.loc[2, protein] = mean_df.get(f"{protein}_dCpt", pd.Series([None])).iloc[0]

df_protein_mean["Value_type"] = ["Tm", "T_opt", "dCpt"]

dump_pickle(df_protein_mean, f"{outdir}/df_protein_mean.pkl")
#print(df_protein_mean)

# for col in std_df.columns:
#     print(f"{col}: {std_df[col].values[0]}")

logging.info("Calculating coefficient of variation")
# Calculate CV directly
cv_values = std_df.values / mean_df.values  # element-wise division

# Create a new DataFrame with the CV values
cv_df = pd.DataFrame(cv_values, columns=std_df.columns).abs()

# Print the coefficient of variation for each parameter
# for col in cv_df.columns:
#     print(f"{col}: {cv_df[col].values[0]}")
dump_pickle(cv_df, f"{outdir}/cv_df.pkl")

cv_df_temp = cv_df.loc[:, cv_df.columns.str.endswith("_Tm") | cv_df.columns.str.endswith("_Topt")]
dump_pickle(cv_df_temp, f"{outdir}/cv_df_temp.pkl")

protein_IDs = pd.read_csv("../data/model_enzyme_params.csv").iloc[:,0].tolist()


highest_r2_row = file.loc[file["r2"].idxmax()]
top_r2_rows = file.nlargest(5, "r2")  # Returns top 5 rows with highest "r2"


best_particle = highest_r2_row.to_frame().T
print(top_r2_rows)


#-------------------------------------------------------------------------------------------------------------------------
logging.info("Retrieving values for T_m and T_opt")

x_values = []
y_values = []  
for ID in protein_IDs:
    tm = ID + "_Tm"
    topt = ID + "_Topt"

    x_values.append(std_df[tm].values[0])
    y_values.append(std_df[topt].values[0])


logging.info("Plotting Standard Deviation (T_m, T_opt) of all proteins")

# Create scatter plot
plt.figure(figsize=(10, 5))
plt.scatter(x_values, y_values, color='b', alpha=0.6)

# Labels and title
plt.xlabel("T_m")
plt.ylabel("T_opt")
plt.title("Scatter Plot Standard Deviation of T_m and T_opt")

# Save and show
plt.savefig("../figures/analysis/scatter_plot_T_m_and_T_opt.png")
plt.show()

#-------------------------------------------------------------------------------------------------------------------------

logging.info("Plotting mean and std for Temperatures")
filtered_columns = [col for col in std_df.columns if not col.endswith("_dCpt")]

x_values = mean_df[filtered_columns].values[0]
y_values = std_df[filtered_columns].values[0]

# Initialize color list
colors = []

# Loop through the columns and assign colors based on the name
for col in filtered_columns:
    if col.endswith("_Topt"):
        colors.append('b')  # Color for Topt
    elif col.endswith("_Tm"):
        colors.append('r')  # Color for Tm
    else:
        colors.append('g')  # Default color for others


# Create scatter plot
plt.figure(figsize=(10, 5))
plt.scatter(x_values, y_values, c =colors, alpha=0.6)

# Labels and title
plt.xlabel("Mean")
plt.ylabel("Std")
plt.title("Scatter Plot Standard Deviation and Mean")

# Get current tick positions
xticks = plt.gca().get_xticks()

# Set new labels (Kelvin to Celsius conversion)
plt.gca().set_xticklabels([f"{t - 273.15:.1f}" for t in xticks])

legend_labels = {
    'Topt': 'b',  # Blue for Topt
    'Tm': 'r',    # Red for Tm
}

handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='b', markersize=10, alpha=0.6, label='Topt'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='r', markersize=10, alpha=0.6, label='Tm')
]
plt.legend(handles=handles)

plt.savefig("../figures/analysis/scatter_plot_mean_std _T.png")

#-------------------------------------------------------------------------------------------------------------------------
logging.info("Plotting mean and std for dCpt")
filtered_columns = [col for col in std_df.columns if col.endswith("_dCpt")]

x_values = mean_df[filtered_columns].values[0]
y_values = std_df[filtered_columns].values[0]



# Create scatter plot
plt.figure(figsize=(10, 5))
plt.scatter(x_values, y_values, color='g', alpha=0.6)

# Labels and title
plt.xlabel("Mean")
plt.ylabel("Std")
plt.title("Scatter Plot Standard Deviation and Mean of dCpt")

plt.savefig("../figures/analysis/scatter_plot_mean_std _dCpt.png")

#-------------------------------------------------------------------------------------------------------------------------
logging.info("Plotting Coefficient of Variation")

x_values = range(cv_df.shape[1])
y_values = cv_df.values

colors = []

for col in cv_df.columns:
    if col.endswith("_Topt"):
        colors.append('b')  # Color for Topt
    elif col.endswith("_Tm"):
        colors.append('r')  # Color for Tm
    else:
        colors.append('g')  # Default color for others


# Create scatter plot
plt.figure(figsize=(10, 5))
plt.scatter(x_values, y_values, c=colors, alpha=0.6)

# Labels and title
plt.ylabel("CV")
plt.title("Scatter Plot Coefficient of Variation")

legend_labels = {
    'Topt': 'b',  # Blue for Topt
    'Tm': 'r',    # Red for Tm
    'dCpt': 'g',  # Green for dCpt
}

handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='b', markersize=10, alpha=0.6, label='Topt'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='r', markersize=10, alpha=0.6, label='Tm'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='g', markersize=10, alpha=0.6, label='dCpt')

]

plt.legend(handles=handles)


plt.savefig("../figures/analysis/scatter_plot_cv.png")

#-------------------------------------------------------------------------------------------------------------------------
logging.info("Plotting mean/mean Tm/Topt")
filtered_columns1 = [col for col in mean_df.columns if col.endswith("_Tm")]
filtered_columns2 = [col for col in mean_df.columns if col.endswith("_Topt")]

x_values = mean_df[filtered_columns1].values[0]
y_values = mean_df[filtered_columns2].values[0]



# Create scatter plot
plt.figure(figsize=(10, 5))
plt.scatter(x_values, y_values, color='g', alpha=0.6)

# Labels and title
plt.xlabel("Mean T_m")
plt.ylabel("Mean T_opt")
plt.title("Scatter Plot Standard Deviation and Mean Tm and Mean Topt")

# Get current tick positions
xticks = plt.gca().get_xticks()
yticks = plt.gca().get_yticks()

# Set new labels (Kelvin to Celsius conversion)
plt.gca().set_xticklabels([f"{t - 273.15:.1f}" for t in xticks])
plt.gca().set_yticklabels([f"{t - 273.15:.1f}" for t in yticks])

plt.savefig("../figures/analysis/scatter_plot_mean_Mean.png")

#-------------------------------------------------------------------------------------------------------------------------
logging.info("Plotting Coefficient of Variation For Temperature only")

x_values = range(cv_df_temp.shape[1])
y_values = cv_df_temp.values

colors = []

for col in cv_df_temp.columns:
    if col.endswith("_Topt"):
        colors.append('b')  # Color for Topt
    elif col.endswith("_Tm"):
        colors.append('r')  # Color for Tm
    else:
        colors.append('g')  # Default color for others


# Create scatter plot
plt.figure(figsize=(10, 5))
plt.scatter(x_values, y_values, c=colors, alpha=0.6)

# Labels and title
plt.ylabel("CV")
plt.title("Scatter Plot Coefficient of Variation")

legend_labels = {
    'Topt': 'b',  # Blue for Topt
    'Tm': 'r',    # Red for Tm
}

handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='b', markersize=10, alpha=0.6, label='Topt'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='r', markersize=10, alpha=0.6, label='Tm')
]

plt.legend(handles=handles)


plt.savefig("../figures/analysis/scatter_plot_cv_temp.png")

logging.info("DONE")