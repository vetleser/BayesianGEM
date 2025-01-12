import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load the UMAP embedding from the pickle file
def load_pickle(filename):
    return pickle.load(open(file=filename, mode='rb'))

outdir = '../results/crowdingDE'
treshold = "R098"



# Load the UMAP embedding from the pickle file

evo_model_frame = load_pickle("../results/crowdingDE/distance_frame.pkl")

#evo_UMAP = load_pickle(f"{outdir}/evo_umap.pkl")
evo_combined_df = load_pickle(f"../results/crowdingDE/evo_combined_df_{treshold}.pkl")

frame_ID = evo_combined_df["frame_ID"].to_numpy()
particle_id = evo_combined_df["particle_ID"].to_numpy()

scaling_factor_values = np.unique(evo_model_frame["scaling_factor"])
crossover_prob_values = np.unique(evo_model_frame["crossover_prob"])
lookup_frame = evo_model_frame.set_index(keys=["scaling_factor","crossover_prob","simulation"])
lookup_frame["frame_ID"] = range(lookup_frame.shape[0])

#All plots
n_neighbors_list = np.array([2, 5, 10, 25, 50, 100, 200])
min_dist_list = np.array([0.1, 0.25])

for n_neighbors in n_neighbors_list:
    for min_dist in min_dist_list:
        embedding = load_pickle(f"../results/crowdingDE/umap/evo_umap_{treshold}_{n_neighbors}_{min_dist}.pkl")

        plt.figure(figsize=(8, 6))


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
            plt.scatter(embedding[this_idxs,0],embedding[this_idxs, 1],edgecolors=["blue","orange", "green", "red", "purple",][simulation],
                        marker=marker,facecolors='none',s=100,label=f"Simulation {simulation+1}")
            plt.title(f"n_neighbors: {n_neighbors}, min_dist: {min_dist}, treshold: {treshold}")
            plt.xlabel("UMAP Dimension 1")
            plt.ylabel("UMAP Dimension 2")
            plt.grid(True)

        plt.legend(loc="upper right")
        plt.savefig(f"../figures/test/umap_evo_{treshold}_{n_neighbors}_{min_dist}.png")

#One plot

# Extract embeddings (assuming umap_ordination is a tuple of (embedding, embedding))
# embedding = evo_UMAP

# Plotting the UMAP results
# plt.figure(figsize=(8, 6))


# marker_dict = {0: 'o', 1: "v", 2: "p",3: "*", 4: "P"} # Based on simulation number
# scaling_factor = 0.5
# crossover_prob = 0.99
# for simulation in range(5):
#     frame_to_choose = lookup_frame.loc[(scaling_factor,crossover_prob,simulation)]
#     this_frame_ID = frame_to_choose["frame_ID"]
#     label_idxs = np.flatnonzero(frame_ID == this_frame_ID)
#     this_idxs = label_idxs
#     marker = marker_dict[simulation]
#     print(label_idxs.shape)
#     plt.scatter(embedding[this_idxs,0],embedding[this_idxs, 1],edgecolors=["blue","orange", "green", "red", "purple",][simulation],
#                 marker=marker,facecolors='none',s=100,label=f"Simulation {simulation+1}")
#     plt.title("UMAP Embedding")
#     plt.xlabel("UMAP Dimension 1")
#     plt.ylabel("UMAP Dimension 2")
#     plt.grid(True)

# plt.legend(loc="upper right")
# plt.savefig("../figures/test/umap_evo.png")

