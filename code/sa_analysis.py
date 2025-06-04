#!/usr/bin/env python
# coding: utf-8

import pickle
import matplotlib
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import multiprocessing
import logging
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')

# Convenient pickle wrappers
def load_pickle(filename):
    return pickle.load(open(file=filename,mode='rb'))

def dump_pickle(obj,filename):
    return pickle.dump(obj=obj,file=open(file=filename, mode='wb'))

def build_a_dataframe_for_all_particles(file, n_priors = 128, r2_threshold = 0.9):
    results = load_pickle(file)

    particles = results.all_particles
    distances = results.all_distances
    
    # Ensure lengths match
    min_len = min(len(particles), len(distances))
    if len(particles) != len(distances):
        logging.warning(f"Length mismatch: {len(particles)} particles vs {len(distances)} distances. Truncating to {min_len}.")
        particles = particles[:min_len]
        distances = distances[:min_len]
    columns = list(results.all_particles[0].keys())
    columns.sort()
    logging.info("Iterating over particles")
    data = list()
    for p in results.all_particles:
        data.append([p[k] for k in columns])
    logging.info("Creating Data Frame")
    df = pd.DataFrame(data=data,columns=columns)
    df['r2'] = results.all_distances
    logging.info(df.shape)
    
    
    logging.info("Doing filtering and labelling of Data Frame")
    df['r2'] = -df['r2']
    df["period"] = "Intermediate"
    df.loc[:n_priors,"period"] = "Prior"
    df.loc[df["r2"] > r2_threshold,"period"] = 'Posterior'
    # Remove samples with a R2 score smaller than -3
    sel_index = df.index[df['r2']>-3]    
    df = df.loc[sel_index,:]
    logging.info(df.shape)

    return df

def inspect_acceptance_rate(filename, n_priors = 128):
    logging.info(f"Loading file: {filename}")
    results = load_pickle(filename)
    acceptance_rates = results.acceptance_rates

    return acceptance_rates 

def count_particles_by_r2(df, filename, combined_count_df):
    r2_col = df['r2']
    
    row = {'file': filename}
    
    for threshold in [0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99]:
        count = (r2_col > threshold).sum()
        logging.info(f"Count of particles with r2 > {threshold}: {count}")
        row[f"r2_{threshold}"] = count

    combined_count_df.loc[len(combined_count_df)] = row  # Append new row
    return count

outdir = '../results/sa'

file_0 = f'{outdir}/smcsa_gem_0.001_0.pkl'
file_1 = "../results/sa/smcsa_gem_may19_0.0001_1.pkl"
file_2 = f'{outdir}/smcsa_gem_may19_0.0001_0.pkl'

file_3 = f'{outdir}/smcsa_gem_may22_0.1_0.pkl'
file_4 = f'{outdir}/smcsa_gem_may22_0.5_0.pkl'
file_5 = f'{outdir}/smcsa_gem_may22_0.5_1.pkl'

file_6 = f'{outdir}/smcsa_gem_may31_0.5_0.pkl'
file_7 = f'{outdir}/smcsa_gem_may31_1.0_0.pkl'
file_8 = f'{outdir}/smcsa_gem_may31_5.0_0.pkl'
file_9 = f'{outdir}/smcsa_gem_may31_10.0_0.pkl'

filenames = ['smcsa_gem_0.001_0', 'smcsa_gem_may19_0.0001_0',
             'smcsa_gem_may22_0.1_0', 'smcsa_gem_may22_0.5_0', 'smcsa_gem_may22_0.5_1',
             'smcsa_gem_may31_0.5_0', 'smcsa_gem_may31_1.0_0', 'smcsa_gem_may31_5.0_0', 'smcsa_gem_may31_10.0_0'] #Removed file_1, different length of all_particles and all_distances due to timeouterror

files = [file_0, file_2, file_3, file_4, file_5, file_6, file_7, file_8, file_9] #Removed file_1, different length of all_particles and all_distances due to timeouterror
# for file, filename in zip(files, filenames):
#     logging.info(f"Processing file: {file}")
#     df = build_a_dataframe_for_all_particles(file)
#     dump_pickle(df, f"{outdir}/{filename}_df.pkl")
#     logging.info(f"Data Frame shape: {df.shape}")
#     logging.info(f"Data Frame columns: {df.columns}")
#     logging.info(f"Data Frame head: {df.head()}")
#     logging.info(f"Data Frame tail: {df.tail()}")
columns = ['file'] + [f'r2_{threshold}' for threshold in [0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99]]
combined_count_df = pd.DataFrame(columns=columns) 

for filename in filenames:
    df = load_pickle(f"{outdir}/{filename}_df.pkl")
    logging.info(f"Loaded Data Frame for {filename} with shape: {df.shape}")
    count_particles_by_r2(df,filename,  combined_count_df)

logging.info("Combined count Data Frame:")
with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 1000):
    logging.warning(f"Combined count Data Frame:\n{combined_count_df}")
filename = '../results/analysis/df_simulation_0_R098.pkl'

logging.info(f"Loading Data Frame from {filename}")
df_0 = load_pickle(filename)
# logging.info(f"Data Frame 0 shape: {df_0.shape}")
# logging.info(f"Data Frame 0 columns: {df_0.columns}")
# logging.info(f"Data Frame 0 head: {df_0.head()}")
# logging.info(f"Data Frame 0 tail: {df_0.tail()}")
#count_particles_by_r2(df_0)
# Uncomment the following lines to build a DataFrame for all particles

# df_0 = build_a_dataframe_for_all_particles(file_0)

# logging.info(f"Data Frame 0 shape: {df_0.shape}")
# logging.info(f"Data Frame 0 columns: {df_0.columns}")
# logging.info(f"Data Frame 0 head: {df_0.head()}")
# logging.info(f"Data Frame 0 tail: {df_0.tail()}")

# acceptance_rates_0 = inspect_acceptance_rate(file_0)
# acceptance_rates_1 = inspect_acceptance_rate(file_1)
# acceptance_rates_2 = inspect_acceptance_rate(file_2)
# acceptance_rates_3 = inspect_acceptance_rate(file_3)
# acceptance_rates_4 = inspect_acceptance_rate(file_4)
# acceptance_rates_5 = inspect_acceptance_rate(file_5)

# plt.figure(figsize=(10, 5))
# plt.plot(acceptance_rates_0, label='Acceptance Rate 0')
# plt.plot(acceptance_rates_1, label='Acceptance Rate 1')
# plt.plot(acceptance_rates_2, label='Acceptance Rate 2')
# plt.plot(acceptance_rates_3, label='Acceptance Rate 3')
# plt.plot(acceptance_rates_4, label='Acceptance Rate 4')
# plt.plot(acceptance_rates_5, label='Acceptance Rate 5')
# plt.xlabel('Generation')
# plt.ylabel('Acceptance Rate')
# plt.title('Acceptance Rate over all Generations')
# plt.legend()
# plt.show()
# plt.savefig('../figures/acceptance_rate.png')
#     df_2 = build_a_dataframe_for_all_particles(file_2)


#--------------------------------------------------------------------

def plot_convergence_inner(distances, maxiter, offset = 128, generation_size = 100, ind_start = 0, label = None):
    # Offset: Number of newly generated particles per generation
    r2s = []
    ind = ind_start
    i = 0
    # This is a cleaver trick to avoid redoing calculations for every iterations which
    # turns out to be very time-consuming.
    r2s_history = -np.array(distances[:offset*maxiter])
    r2s_history_argsorted = np.argsort(r2s_history)
    while ind < len(distances):
        i += 1
        if i > maxiter:
            break
        # This is a mask ensuring data created after the interation are excluded
        filter_mask = r2s_history_argsorted < ind + offset
        r2s_now = r2s_history[r2s_history_argsorted[filter_mask]][-generation_size:]
        r2s.append(r2s_now)
        ind += offset
    y = np.array(r2s)
    yp = np.percentile(y,[5,50,95],axis=1)
    plt.plot(np.arange(len(r2s)),yp[1,:], label = label)
    plt.fill_between(np.arange(len(r2s)),yp[0,:],yp[2,:],alpha=0.5)

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 30}


# model_frame = load_pickle("../results/sa/distance_frame.pkl")
# logging.info(f"Model frame shape: {model_frame.shape}")
# logging.info(f"Model frame columns: {model_frame.columns}")
# logging.info(f"Model frame head: {model_frame.head()}")

# for idx, row in model_frame.iterrows():
#     logging.info(f"Processing row index: {idx}")
#     # logging.info(f"Processing row: {row.Index}")
#     # logging.info(f"Final temperature: {row.final_temp}, Move type: {row.move_type}, Step size: {row.step_size}, Simulation: {row.simulation}")

#     #row = model_frame.iloc[0]
#     all_distances = row["all_distances"]
#     population = row["population"]
#     indices = [i for i, p in enumerate(population[-1]) if all_distances[p] < -0.97]
#     final_distances = [all_distances[p] for p in population[-1] if all_distances[p] < 3]
#     logging.info(f"Final generation distances: {final_distances}")
#     mean_distances = []
#     gen_iter = 0
#     for gen in population:
#         gen_iter += 1
#         #logging.info(f"Processing generation {gen_iter} with {len(gen)} particles")
#         # distances = [all_distances[p] for i, p in enumerate(gen) if all_distances[p] < 5 and i in indices]
#         distances = [all_distances[p] for i, p in enumerate(gen) if i in indices]

#         logging.info(f"Generation {gen_iter} n_distances: {len(distances)}")
#         mean_distances.append(-np.mean(distances))

#     plt.figure(figsize=(10, 5))
#     plt.plot(mean_distances, label='Mean Distances')
#     #plt.ylim([-2, 1])
#     plt.xlabel('Generation')
#     plt.ylabel('Mean Distance')
#     plt.title(f'Mean Distances over Generations Simulation {idx}')
#     plt.legend()
#     plt.savefig(f'../figures/mean_distances_{idx}.png')

# matplotlib.rc('font', **font)
# proper_names = {'unpermuted': "Unpermuted", 'permuted_0': "Permuted 1",
#                 'permuted_1': "Permuted 2", 'permuted_2': "Permuted 3"}
# maxiter = 1000
# i = 1
# plt.figure(figsize=(20,20))
# for index, series in model_frame["all_distances"]:
#     plt.subplot(2,2,i)
#     origin_distances = series[index]
#     # Simulation 1
#     original_distances = origin_distances["original"]
#     # Simulation 2
#     replicate_distances = origin_distances["replicate"]
#     plot_convergence_inner(original_distances,maxiter, label = 'Simulation 1')
#     plot_convergence_inner(replicate_distances,maxiter, label = 'Simulation 2')
#     if i==4:
#         handles, labels = plt.gca().get_legend_handles_labels()
#     plt.ylim([0,1])
#     plt.xlabel('Iterations')
#     plt.ylabel('$R^2$')
#     plt.title(proper_names[index])
#     i += 1
#     plt.tight_layout()
# plt.subplots_adjust(bottom=0.1)
# plt.gcf().legend(handles,labels, loc=(.34,0.005),ncol=2,handletextpad=0.5)
# plt.savefig("../figures/R2_vetle.pdf")
# plt.show()
logging.info("DONE")

