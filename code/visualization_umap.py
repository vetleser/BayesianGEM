import pickle
import matplotlib.pyplot as plt
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

#Individual plots
# n_neighbors_list = np.array([2, 5, 10, 25, 50, 100, 200])
# min_dist_list = np.array([0.0, 0.05])

# for n_neighbors in n_neighbors_list:
#     for min_dist in min_dist_list:
#         embedding = load_pickle(f"../results/crowdingDE/umap/evo_umap_{treshold}_{n_neighbors}_{min_dist}.pkl")

#         plt.figure(figsize=(8, 6))


#         marker_dict = {0: 'o', 1: "v", 2: "p",3: "*", 4: "P"} # Based on simulation number
#         scaling_factor = 0.5
#         crossover_prob = 0.99
#         for simulation in range(5):
#             frame_to_choose = lookup_frame.loc[(scaling_factor,crossover_prob,simulation)]
#             this_frame_ID = frame_to_choose["frame_ID"]
#             label_idxs = np.flatnonzero(frame_ID == this_frame_ID)
#             this_idxs = label_idxs
#             marker = marker_dict[simulation]
#             print(label_idxs.shape)
#             plt.scatter(embedding[this_idxs,0],embedding[this_idxs, 1],edgecolors=["blue","orange", "green", "red", "purple",][simulation],
#                         marker=marker,facecolors='none',s=100,label=f"Simulation {simulation+1}")
#             plt.title(f"n_neighbors: {n_neighbors}, min_dist: {min_dist}, treshold: {treshold}")
#             plt.xlabel("UMAP Dimension 1")
#             plt.ylabel("UMAP Dimension 2")
#             plt.grid(True)

#         plt.legend(loc="upper right")
#         plt.savefig(f"../figures/test/umap_evo_{treshold}_{n_neighbors}_{min_dist}.png")

#One plot
n_neighbors_list = np.array([10, 25, 50, 100, 200])
min_dist_list = np.array([0.0, 0.05, 0.1, 0.25])

n_rows = len(n_neighbors_list)
n_columns = len(min_dist_list)

#plt.figure(figsize=(8*n_columns, 6*n_rows))
fig, axes = plt.subplots(n_rows, n_columns, figsize=(8*n_columns, 6*n_rows), constrained_layout=False)

plt.subplots_adjust(top=0.95, left=0.1, right=0.95, bottom=0.05, wspace=0, hspace=0)

marker_dict = {0: 'o', 1: "v", 2: "p",3: "*", 4: "P"} # Based on simulation number
scaling_factor = 0.5
crossover_prob = 0.99

# fig.text(8*n_columns/2, 10, "min_dist", rotation='horizontal', fontsize=100)
# fig.text(10, 6*n_rows/2, "n_neighbors", rotation='vertical', fontsize=100)

fig.text(0.5, 1-0.02, "min_dist", ha='center', fontsize=48)  # X-axis at bottom
fig.text(0.02, 0.5, "n_neighbors", va='center', rotation='vertical', fontsize=48)  # Y-axis at left

fig.text(0.05, 0.95, "A" if treshold=="R098" else "B", ha='center', fontsize=98, fontweight="bold")

for row, n_neighbors in enumerate(n_neighbors_list):
    for col, min_dist in enumerate(min_dist_list):
        embedding = load_pickle(f"../results/crowdingDE/umap/evo_umap_{treshold}_{n_neighbors}_{min_dist}.pkl")

        ax = axes[row, col] if n_rows > 1 and n_columns > 1 else axes[row if n_rows > 1 else col]

        for simulation in range(5):
            frame_to_choose = lookup_frame.loc[(scaling_factor,crossover_prob,simulation)]
            this_frame_ID = frame_to_choose["frame_ID"]
            label_idxs = np.flatnonzero(frame_ID == this_frame_ID)
            this_idxs = label_idxs
            marker = marker_dict[simulation]
            print(label_idxs.shape)
            ax.scatter(embedding[this_idxs,0],embedding[this_idxs, 1],edgecolors=["blue","orange", "green", "red", "purple",][simulation],
                        marker=marker,facecolors='none',s=100,label=f"Simulation {simulation+1}")
            #plt.title(f"n_neighbors: {n_neighbors}, min_dist: {min_dist}, treshold: {treshold}")
            ax.grid(True)
        
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(2)  # Change thickness to 2 (or any value you like)


        if(row == 0):
            ax.set_title(min_dist, fontsize= 36)

        if(col == 0):
            ax.set_ylabel(n_neighbors, fontsize = 36)
            

print("DONE")
#plt.savefig(f"../figures/evo_umap_all_{treshold}_withletter.png")



