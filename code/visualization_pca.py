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

evo_model_frame = load_pickle("../results/crowdingDE/distance_frame.pkl")

def plot_convergence_evo(distances, maxiter, populations, label = None):
    distance_array = -np.array(distances)
    yp = np.vstack([np.percentile(distance_array[population],[5,50,95]) for population in populations[0:maxiter]])
    plt.plot(np.arange(yp.shape[0]),yp[:,1], label = label)
    plt.fill_between(np.arange(yp.shape[0]),yp[:,0],yp[:,2],alpha=0.5)

evo_PCS,evo_EV = load_pickle("../results/crowdingDE/evo_pca.pkl")

evo_combined_df = load_pickle("../results/crowdingDE/evo_combined_df.pkl")

frame_ID = evo_combined_df["frame_ID"].to_numpy()
particle_id = evo_combined_df["particle_ID"].to_numpy()

scaling_factor_values = np.unique(evo_model_frame["scaling_factor"])
crossover_prob_values = np.unique(evo_model_frame["crossover_prob"])
lookup_frame = evo_model_frame.set_index(keys=["scaling_factor","crossover_prob","simulation"])
lookup_frame["frame_ID"] = range(lookup_frame.shape[0])

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 30}

matplotlib.rc('font', **font)
maxiter = 1000
plt.figure(figsize=(24,12))
plt.subplot(1,2,1)
for _, entry in evo_model_frame.iterrows():
    scaling_factor = entry["scaling_factor"]
    crossover_prob = entry["crossover_prob"]
    if scaling_factor != 0.5 or crossover_prob != 0.99:
        continue
    simulation = entry["simulation"]
    populations = entry["population"]
    distances = entry["all_distances"]
    plot_convergence_evo(distances,maxiter=1000,populations=populations, label = f'Simulation {simulation + 1}')
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.ylim([0,1])
    plt.xlabel('Generation')
    plt.ylabel('$R^2$')
    plt.title("A",loc="left",fontsize=45,fontweight="bold")
    plt.tight_layout()
    plt.legend(loc="lower right")



plt.subplot(1,2,2)
marker_dict = {0: 'o', 1: "v", 2: "p",3: "*", 4: "P"} # Based on simulation number
scaling_factor = 0.5
crossover_prob = 0.99
for simulation in range(5):
    frame_to_choose = lookup_frame.loc[(scaling_factor,crossover_prob,simulation)]
    this_frame_ID = frame_to_choose["frame_ID"]
    label_idxs = np.flatnonzero(frame_ID == this_frame_ID)
    this_idxs = label_idxs
    marker = marker_dict[simulation]
    print(label_idxs.shape)
    plt.scatter(evo_PCS[this_idxs,0],evo_PCS[this_idxs, 1],edgecolors=["blue","orange", "green", "red", "purple",][simulation],
                marker=marker,facecolors='none',s=100,label=f"Simulation {simulation+1}")
    plt.xlabel('PC1 ({:.2f}%)'.format(evo_EV[0]*100))
    plt.ylabel('PC2 ({:.2f}%)'.format(evo_EV[1]*100))
    plt.xticks(np.arange(-15,20,5))
    plt.yticks(np.arange(-20,60,10))
    plt.xlim((-17,10))
    plt.ylim((-24,40))
    plt.title("B",loc="left",fontsize=45,fontweight="bold")
    plt.tight_layout()
    colors = ['#1f78b4', 'grey', '#fc8d59']

plt.legend(loc="upper right")
plt.savefig("../figures/evo_fig.png",dpi=400)
plt.show()