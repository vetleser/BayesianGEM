import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.cluster as cluster
import pickle
import GEMS
import itertools
import dill
import abc_etc
import evo_etc

# Convenient pickle wrappers
def load_pickle(filename):
    return dill.load(open(file=filename,mode='rb'))
        
def dump_pickle(obj,filename):
    return dill.dump(obj=obj,file=open(file=filename, mode='wb'))

'r_0961No1' # Pyruvate dehydrogenase
['r_0959No1','r_0959No2','r_0959No3'] # Pyruvate decarboxylase
'r_0450No1' # Fructose-bisphosphate aldolase
'r_0438No1'# Ferrocytochrome-c:oxygen oxidoreductase
'r_0917No1' # Phosphoserine phosphatase
'r_0997No1' # Shikimate kinase
'r_2111' # Growth



signature_reactions = {'PDH': 'r_0961No1', 'FBA': 'r_0450No1', 'FCO': 'r_0438No1', 'PSP': 'r_0917No1', 'SHK': 'r_0997No1', 'GRW': 'r_2111'}
signature_full_name = {'PDH': 'Pyruvate dehydrogenase', 'FBA': 'Fructose-bisphosphate aldolase', 'FCO': 'Ferrocytochrome-c:oxygen oxidoreductase',
                       'PSP': 'Phosphoserine phosphatase', 'SHK': 'Shikimate kinase', 'GRW': "Growth"}

fva_res = load_pickle("../results/crowdingDE/evo_fva_R098.pkl")

def aggregate_fva_results(result_df,simulation_attributes):
    flattened_df_list = []
    for _, row in result_df.drop(columns=["particle"]).iterrows():
        raw_df = row["fva_res"]
        df = raw_df
        for attribute in simulation_attributes:
            df[attribute] = row[attribute]
        flattened_df_list.append(df)

    combined_fva_frame = (
        pd.concat(flattened_df_list).
        assign(range= lambda df: df["maximum"] - df["minimum"],
                                                            midpoint= lambda df: (df["maximum"] + df["minimum"]) / 2).
        drop(columns=["minimum", "maximum"])
            )
    simulation_attributes.extend(["condition","reaction","T"])
    aggregated_fva_res = (
        combined_fva_frame.replace([np.inf, -np.inf],np.nan).
        dropna(how="all").
        groupby(simulation_attributes).
        agg(["mean","min","max","std","count"])
                        )
    return aggregated_fva_res

aggregated_fva_res = aggregate_fva_results(fva_res, ["simulation"])

aggregated_fva_res.reset_index(inplace=True)

aggregated_fva_res.reset_index()

reactions = signature_reactions.keys()
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 45}
matplotlib.rc('font', **font)
T_0 = 273.15
def extract_results_and_plot(simulation,reaction, what, linestyle,condition="aerobic",label=None):
    if label is None:
        label = f"Simulation {simulation+1}"
    react_id = signature_reactions[reaction]
    subsetted_frame = aggregated_fva_res.reset_index().pipe(lambda df: df[(df["simulation"] == simulation) & (df["reaction"] == react_id) & (df["condition"] == condition)])
    T = subsetted_frame[("T","")]
    mu = subsetted_frame[(what,"mean")]
    minimum = subsetted_frame[(what,"min")]
    maximum = subsetted_frame[(what,"max")]
    plt.errorbar(T-T_0,mu,yerr=np.row_stack((mu-minimum,maximum-mu)),markersize=8,capsize=15,linestyle=linestyle, label=label)
    return
    
plt.figure(figsize=(63,60))
nrows = 4
ncols = 3
i = 1
subplot_order = {1: 1, 2: 4, 3: 2, 4: 5, 5: 3, 6: 6, 7: 7, 8: 10, 9: 8, 10: 11, 11: 9, 12: 12}
ymaxs = {"PDH": 1, "FBA": 30, "FCO": 10, "PSP": 0.5, "SHK": 0.2, 'GRW': 0.4 }
for reaction in reactions:
    for what in ["midpoint", "range"]:
        if i == 12:
            # The growth range is not interesting, so we use it for legend instead
            continue
        plt.subplot(4,3,subplot_order[i])
        extract_results_and_plot(simulation=0,reaction=reaction, what=what, linestyle="solid")
        extract_results_and_plot(simulation=1,reaction=reaction, what=what, linestyle="solid")
        extract_results_and_plot(simulation=2,reaction=reaction, what=what, linestyle="solid")
        extract_results_and_plot(simulation=3,reaction=reaction, what=what, linestyle="solid")
        extract_results_and_plot(simulation=4,reaction=reaction, what=what, linestyle="solid")
        plt.ylabel(r"Flux range $\left[\frac{\mathrm{mmol}}{\mathrm{gDW}\cdot \mathrm{h}}\right]$" if what == "range"
               else r"Flux midpoint $\left[\frac{\mathrm{mmol}}{\mathrm{gDW}\cdot \mathrm{h}}\right]$")
        plt.xlabel(r"Temperature $\left[^\circ\mathrm{C}\right]$")
        plt.title(f"{signature_full_name[reaction]} {what}")
        plt.xlim((15,45))
        plt.ylim((-0.01,ymaxs[reaction]))
        if i == 10:
            handles, labels = plt.gca().get_legend_handles_labels()
        i += 1
        plt.tight_layout()
plt.gcf().legend(handles,labels, loc=(.7,0.1),ncol=1,handletextpad=0.5)
plt.savefig("../figures/evo_aerobic_fva_R098_vetle.png")