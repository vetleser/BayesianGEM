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


logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')


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
    


outdir = "../results/toy_example"
n_simulations = 4
n_plots = 4

toy_example_results = pd.DataFrame(index=pd.MultiIndex.from_product([range(n_plots),range(n_simulations)],names=["Plot","Simulation"])).reset_index()
toy_example_results["modelfile"] = list(itertools.starmap(lambda plot,simulation: f"{'simanneal'}_{plot}_{simulation}.pkl",
                                                       toy_example_results[["Plot","Simulation"]].itertuples(index=False,name=None)))
toy_example_results["model"] = list(map(lambda filename: load_pickle(f"{outdir}/{filename}"),toy_example_results["modelfile"]))
toy_example_results["final_generation"] = list(itertools.starmap(lambda method, model: ( extract_simulated_annealing_final_generation)(model),
                                               toy_example_results[["Plot","model"]].itertuples(index=False,name=None)))
toy_example_results.set_index(["Plot","Simulation"],inplace=True)

logging.info("Toy example results loaded")

particles = toy_example_results.loc[(3, 3), "final_generation"]

logging.info(particles)
tol = 0.5

def sort_particles(particles, plot, sim):
    close_to_minima = [0, 0, 0, 0]
    minima = [(-1, 1), (1, 1), (-1, -1), (1, -1)]
    for p in particles:
        for i, m in enumerate(minima):
            if (p["x"] - m[0])**2 + (p["y"] - m[1])**2 < tol**2:
                close_to_minima[i] += 1

    logging.info(f"Close to minima in plot {plot}, sim {sim}: {close_to_minima}")
    logging.info(f"Total close to minima: {sum(close_to_minima)}")

#sort_particles(particles)

for j in range(n_plots):
    for i in range(n_simulations):
        particles = toy_example_results.loc[(j, i), "final_generation"]
        sort_particles(particles, j, i)
        #logging.info(f"Plot {j}, Simulation {i}: {particles}")
        # dump_pickle(particles, f"{outdir}/particles_{i}_{j}.pkl")
        # dump_pickle(toy_example_results.loc[(i, j), "model"], f"{outdir}/model_{i}_{j}.pkl")

logging.info("Various minima results")
toy_example_results = pd.DataFrame(index=pd.MultiIndex.from_product([range(n_plots),range(n_simulations)],names=["Plot","Simulation"])).reset_index()
toy_example_results["modelfile"] = list(itertools.starmap(lambda plot,simulation: f"{'simanneal_variousminima'}_{plot}_{simulation}.pkl",
                                                       toy_example_results[["Plot","Simulation"]].itertuples(index=False,name=None)))
toy_example_results["model"] = list(map(lambda filename: load_pickle(f"{outdir}/{filename}"),toy_example_results["modelfile"]))
toy_example_results["final_generation"] = list(itertools.starmap(lambda method, model: ( extract_simulated_annealing_final_generation)(model),
                                               toy_example_results[["Plot","model"]].itertuples(index=False,name=None)))
toy_example_results.set_index(["Plot","Simulation"],inplace=True)


particles = toy_example_results.loc[(3, 3), "final_generation"]

for j in range(n_plots):
    for i in range(n_simulations):
        particles = toy_example_results.loc[(j, i), "final_generation"]
        sort_particles(particles, j, i)
        #logging.info(f"Plot {j}, Simulation {i}: {particles}")
        # dump_pickle(particles, f"{outdir}/particles_{i}_{j}.pkl")
        # dump_pickle(toy_example_results.loc[(i, j), "model"], f"{outdir}/model_{i}_{j}.pkl")

