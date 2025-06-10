#!/usr/bin/env python
# coding: utf-8


import itertools
import random_sampler
import numpy as np
import sa_etc 
import os
import math
import logging
import copy
import time
import dill
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')

logging.info("BEGIN")


# Convenient pickle wrappers
def load_pickle(filename):
    return dill.load(open(file=filename,mode='rb'))
        
def dump_pickle(obj,filename):
    return dill.dump(obj=obj,file=open(file=filename, mode='wb'))

def extract_simulated_annealing_final_generation(model: sa_etc.SimulatedAnnealing):
    if not hasattr(model, "population"):
        logging.info("No population in model")
        return []
    if not hasattr(model, "all_particles"):
        logging.info("No all_particles in model")
        return []
    if not model.population:
        logging.info("Empty population")
        return []
    particles_to_pick = model.population[-1]
    return [model.all_particles[particle] for particle in particles_to_pick]

def extract_sa_population(model: sa_etc.SimulatedAnnealing):
    if not hasattr(model, "population"):
        logging.info("No population in model")
        return []
    if not model.population:
        logging.info("Empty population")
        return []
    return model.population  # Return the last population
    
def extract_final_distances(model: sa_etc.SimulatedAnnealing):
    if not hasattr(model, "all_particles"):
        logging.info("No all_particles in model")
        return []
    if not hasattr(model, "all_distances"):
        logging.info("No all_distances in model")
        return []
    particles_to_pick = model.population[-1]
    return [model.all_distances[particle] for particle in particles_to_pick]

def extract_all_distances(model: sa_etc.SimulatedAnnealing):
    if not hasattr(model, "all_particles"):
        logging.info("No all_particles in model")
        return []
    if not hasattr(model, "all_distances"):
        logging.info("No all_distances in model")
        return []
    return [model.all_distances]



outdir = "../results/toy_example"
n_simulations = 4
n_plots = 1
model_types = ["rastr", "gaussian"]

# Create all combinations
index_tuples = list(itertools.product(model_types, range(n_plots), range(n_simulations)))
toy_example_results = pd.DataFrame(index_tuples, columns=["ModelType", "Plot", "Simulation"])

# Construct modelfile
toy_example_results["modelfile"] = toy_example_results.apply(
    lambda row: f"simanneal_{row.ModelType}_{row.Plot}_{row.Simulation}.pkl", axis=1
)

# Load models, populations, distances
toy_example_results["model"] = toy_example_results["modelfile"].map(lambda filename: load_pickle(f"{outdir}/{filename}"))
toy_example_results["population"] = toy_example_results["model"].map(extract_sa_population)
toy_example_results["all_distances"] = toy_example_results["model"].map(extract_all_distances)

# toy_example_results = pd.DataFrame(index=pd.MultiIndex.from_product([range(n_plots),range(n_simulations)],names=["Plot","Simulation"])).reset_index()
# toy_example_results["modelfile"] = list(itertools.starmap(lambda plot,simulation: f"{'simanneal'}_rastr_{plot}_{simulation}.pkl",
#                                                        toy_example_results[["Plot","Simulation"]].itertuples(index=False,name=None)))
# toy_example_results["model"] = list(map(lambda filename: load_pickle(f"{outdir}/{filename}"),toy_example_results["modelfile"]))
# toy_example_results["population"] = list(itertools.starmap(lambda method, model: ( extract_sa_population)(model),
#                                                toy_example_results[["Plot","model"]].itertuples(index=False,name=None)))
# toy_example_results["all_distances"] = list(itertools.starmap(lambda method, model: (extract_all_distances)(model),
#                                                toy_example_results[["Plot","model"]].itertuples(index=False,name=None)))
# toy_example_results.set_index(["Plot","Simulation"],inplace=True)



# logging.info("Toy example results loaded")
with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # Set display options to show all rows and columns
    logging.info("Toy example results:")
    logging.info(toy_example_results)



tol = 0.3

toy_example_setup = load_pickle(f"{outdir}/toy_example_df.pkl")

def sort_particles_by_xy(particles, plot, sim):
    logging.info(f"Number of particles in plot {plot}, simulation {sim}: {len(particles)}")
    minimas = toy_example_setup["Minima_list"].to_list()
    close_to_minima = []
    for particle in particles:
        x = particle["x"]
        y = particle["y"]
        close = False
        for minima in minimas:
            if math.sqrt((x - minima[0])**2 + (y - minima[1])**2) < tol:
                scale = toy_example_setup["Scaling_list"].to_list()[minimas.index(minima)]
                logging.info(f"Particle {particle} is close to minima {minima} with scale  {scale} in plot {plot}, simulation {sim}")
                break
        close_to_minima.append(close)

def sort_particles_by_r2(particles, distances, plot, sim):
    logging.info(f"Number of particles in plot {plot}, simulation {sim}: {len(particles)}")
    # for i, particle in enumerate(particles):
    #     d = distances[i]
    #     logging.info(f"Particle {particle} distance: {d}")

    tresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    for t in tresholds:
        interval = t + 0.1
        counter = 0
        for i, particle in enumerate(particles):
            d = -distances[i]
            if d > t and d < interval:
                counter += 1
        logging.info(f"Plot {plot}, Simulation {sim}, Treshold {t}: {counter} particles with distance  {t} < d < {t + 0.1}")

    # logging.info(f"Close to minima in plot {plot}, sim {sim}: {close_to_minima}")
    # logging.info(f"Total close to minima: {sum(close_to_minima)}")

#sort_particles(particles)

# for j in range(n_plots):
#     for i in range(n_simulations):
#         particles = toy_example_results.loc[(j, i), "final_generation"]
#         distances = toy_example_results.loc[(j, i), "final_distances"]
#         sort_particles_by_r2(particles, distances, j, i)
#         sort_particles_by_xy(particles, j, i)
#         #logging.info(f"Plot {j}, Simulation {i}: {particles}")
#         # dump_pickle(particles, f"{outdir}/particles_{i}_{j}.pkl")
#         # dump_pickle(toy_example_results.loc[(i, j), "model"], f"{outdir}/model_{i}_{j}.pkl")

def plot_convergence(distances, maxiter, label = None):
    distance_array = -np.array(distances)  # shape: (n_generations, n_particles)
    distance_array = distance_array[:maxiter]
    yp = np.percentile(distance_array, [5, 50, 95], axis=1)  # axis=1 = per generation
    plt.plot(np.arange(yp.shape[1]), yp[1], label=label)  # 50th percentile
    plt.fill_between(np.arange(yp.shape[1]), yp[0], yp[2], alpha=0.5)  # 5–95%




font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 30}
matplotlib.rc('font', **font)


logging.info(f"Index: {toy_example_results.index}")

for model_type, group in toy_example_results.groupby("ModelType"):
    plt.figure(figsize=(20, 20))
    if model_type == "rastr":
        plt.suptitle("Convergence of Rastrigin Toy Problem", fontsize=30, fontweight='bold')
    elif model_type == "gaussian":
        plt.suptitle("Convergence of Gaussian Toy Problem", fontsize=30, fontweight='bold')
    
    for _, entry in group.iterrows():
        simulation = entry["Simulation"]
        populations = entry["population"]
        all_distances = entry["all_distances"][0]

        distances = []
        for population in populations:
            distances.append([all_distances[particle] for particle in population])

        plt.subplot(2, 2, simulation + 1)
        plot_convergence(distances, maxiter=200, label=f'Simulation {simulation + 1}')
        plt.ylim([0, 1])
        plt.xlabel('Generation')
        plt.ylabel('$R^2$')
        plt.title(f"Simulation {simulation}")
        plt.tight_layout()
    

    plt.subplots_adjust(bottom=0.10)
    plt.savefig(f"../figures/toy_example/sa_{model_type}_R2.png", dpi=300)
    plt.show()

for model_type, group in toy_example_results.groupby("ModelType"):
    plt.figure(figsize=(20, 20))
    if model_type == "rastr":
        plt.suptitle("Distribution of Distances in Final Population — Rastrigin Toy Example", fontsize=30, fontweight='bold')
    elif model_type == "gaussian":
        plt.suptitle("Distribution of Distances in Final Population — Gaussian Toy Example", fontsize=30, fontweight='bold')

    for _, entry in group.iterrows():
        simulation = entry["Simulation"]
        final_population = entry["population"][-1]
        all_distances = entry["all_distances"][0]

        distances = []
        distances.append([-all_distances[particle] for particle in final_population])

        plt.subplot(2, 2, simulation + 1)
        plt.hist(distances, bins=20, alpha=0.7)
        plt.xlim([0, 1])
        plt.xlabel('$R^2$')
        plt.ylabel('Number of particles')
        plt.title(f"Simulation {simulation}")
        plt.tight_layout()
    plt.subplots_adjust(bottom=0.10)
    plt.savefig(f"../figures/toy_example/sa_{model_type}_histogram.png", dpi=300)
    plt.show()

for model_type, group in toy_example_results.groupby("ModelType"):
    plt.figure(figsize=(24, 16))
    if model_type == "rastr":
        plt.title("Distribution of Distances in Final Population — Rastrigin Toy Example", fontsize=20, fontweight='bold')
    elif model_type == "gaussian":
        plt.title("Distribution of Distances in Final Population — Gaussian Toy Example", fontsize=20, fontweight='bold')

    bins = np.linspace(0, 1, 21)  # 20 bins between 0 and 1

    all_sim_distances = []
    labels = []
    for _, entry in group.iterrows():
        simulation = entry["Simulation"]
        final_population = entry["population"][-1]
        all_distances = entry["all_distances"][0]

        distances = [-all_distances[particle] for particle in final_population]
        all_sim_distances.append(distances)
        labels.append(f'Simulation {simulation}')

    plt.hist(all_sim_distances, bins=20, stacked=True, label=labels)
    plt.xlim([0, 1])
    plt.xlabel('$R^2$')
    plt.ylabel('Number of particles')
    plt.legend(title='Simulation', fontsize='small', title_fontsize='medium', handlelength=1, handleheight=0.7)
    plt.tight_layout()
    plt.savefig(f"../figures/toy_example/sa_{model_type}_combined_histogram.png", dpi=300)
    plt.show()
    
logging.info("DONE")

