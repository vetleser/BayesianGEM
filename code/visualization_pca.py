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

evo_PCS1,evo_EV1 = load_pickle("../results/crowdingDE/evo_pca_R098.pkl")

treshold1 = "R098"
evo_combined_df_R098 = load_pickle(f"../results/crowdingDE/evo_combined_df_{treshold1}.pkl")

frame_ID1 = evo_combined_df_R098["frame_ID"].to_numpy()
particle_id1 = evo_combined_df_R098["particle_ID"].to_numpy()

evo_PCS2,evo_EV2 = load_pickle("../results/crowdingDE/evo_pca_R0985.pkl")


treshold2 = "R0985"
evo_combined_df_R0985 = load_pickle(f"../results/crowdingDE/evo_combined_df_{treshold2}.pkl")

frame_ID2 = evo_combined_df_R0985["frame_ID"].to_numpy()
particle_id2 = evo_combined_df_R0985["particle_ID"].to_numpy()

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
marker_dict = {0: 'o', 1: "v", 2: "p",3: "*", 4: "P"} # Based on simulation number
scaling_factor = 0.5
crossover_prob = 0.99
for simulation in range(5):
    frame_to_choose = lookup_frame.loc[(scaling_factor,crossover_prob,simulation)]
    this_frame_ID = frame_to_choose["frame_ID"]
    label_idxs = np.flatnonzero(frame_ID1 == this_frame_ID)
    this_idxs = label_idxs
    marker = marker_dict[simulation]
    print(label_idxs.shape)
    plt.scatter(evo_PCS1[this_idxs,0],evo_PCS1[this_idxs, 1],edgecolors=["blue","orange", "green", "red", "purple",][simulation],
                marker=marker,facecolors='none',s=100,label=f"Simulation {simulation+1}")
    plt.xlabel('PC1 ({:.2f}%)'.format(evo_EV1[0]*100))
    plt.ylabel('PC2 ({:.2f}%)'.format(evo_EV1[1]*100))
    plt.xticks(np.arange(-20,50,10))
    plt.yticks(np.arange(-20,60,10))
    # plt.xlim((-17,10))
    # plt.ylim((-24,40))
    plt.xlim((-25,44))
    plt.ylim((-19,43))
    plt.text(0, 1.08, "A", transform=plt.gca().transAxes, fontsize=45, fontweight="bold", va="top", ha="left")
    plt.title("R² > 0.98", fontsize=28, fontweight="bold")    
    plt.tight_layout()
    colors = ['#1f78b4', 'grey', '#fc8d59']
    plt.legend(loc="upper right")


plt.subplot(1,2,2)
marker_dict = {0: 'o', 1: "v", 2: "p",3: "*", 4: "P"} # Based on simulation number
scaling_factor = 0.5
crossover_prob = 0.99
for simulation in range(5):
    frame_to_choose = lookup_frame.loc[(scaling_factor,crossover_prob,simulation)]
    this_frame_ID = frame_to_choose["frame_ID"]
    label_idxs = np.flatnonzero(frame_ID2 == this_frame_ID)
    this_idxs = label_idxs
    marker = marker_dict[simulation]
    print(label_idxs.shape)
    plt.scatter(evo_PCS2[this_idxs,0],evo_PCS2[this_idxs, 1],edgecolors=["blue","orange", "green", "red", "purple",][simulation],
                marker=marker,facecolors='none',s=100,label=f"Simulation {simulation+1}")
    plt.xlabel('PC1 ({:.2f}%)'.format(evo_EV2[0]*100))
    plt.ylabel('PC2 ({:.2f}%)'.format(evo_EV2[1]*100))
    plt.xticks(np.arange(-20,50,10))
    plt.yticks(np.arange(-20,60,10))
    # plt.xlim((-17,10))
    # plt.ylim((-24,40))
    plt.xlim((-25,44))
    plt.ylim((-19,43))
    plt.text(0, 1.08, "B", transform=plt.gca().transAxes, fontsize=45, fontweight="bold", va="top", ha="left")
    plt.title("R² > 0.985", fontsize=28, fontweight="bold")
    plt.tight_layout()
    colors = ['#1f78b4', 'grey', '#fc8d59']
    plt.legend(loc="upper right")


plt.savefig(f"../figures/evo_pca_both_fig_new.png",dpi=400)
plt.show()