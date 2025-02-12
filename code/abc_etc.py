#!/usr/bin/env python
# coding: utf-8

# #### Add constraint: Topt<Tm

# In[1]:


import multiprocessing
import numpy as np
import scipy.stats as ss
from multiprocessing import Process,cpu_count,Manager
import pickle
from random_sampler import RV
import logging
from typing import Callable, Dict, Iterable, List
import numpy.typing as npt


simResultType = Dict[str, Dict[str, npt.NDArray[np.float64]]]
candidateType = Dict[str, float]
priorType = Dict[str, RV]
distanceArgType = Dict[str, npt.NDArray[np.float64]]




# In[ ]:

class SMCABC:
    def __init__(self,simulator: Callable[[candidateType], simResultType],priors: priorType,min_epsilon: float,population_size: int,
    distance_function: Callable[[distanceArgType, distanceArgType], float],
                 Yobs: distanceArgType,outfile: str,cores: int=cpu_count(),generation_size: int =128, maxiter: int=100000):
        '''
        simulator:       a function that takes a dictionary of parameters as input. Ouput {'data':Ysim}
        priors:          a dictionary which use id of parameters as keys and RV class object as values
        min_epsilon:     minimal epsilon
        population_size: the size of each population
        distance_function: a function that calculate the distance between observed data and simulated data
        Yobs:            observed data
        outfile:         unique id for the experiment. This will be also used to continue a simulation that 
                         is partly done
        cores:           number of treads
        maxiter:         Maximum number of iterations before breaking the fitting
        
        !!!Important: distance is to be minimized!!!
        '''
        self.simulator = simulator
        self.priors = priors
        self.posterior = priors                   # to be updated
        self.population_size = population_size
        self.distance_function = distance_function
        self.min_epsilon = min_epsilon
        self.Yobs = Yobs
        self.outfile = outfile
        self.population: List[candidateType] = []  # a list of the current population
        self.distances = []    # a list of distances for particles in population
        self.simulations = 0  # number of simulations performed 
        self.cores = cores    
        self.simulated_data_t0 = [] # first population
        self.population_t0 = []
        self.distances_t0 = []
        self.simulated_data = []  # last population
        self.epsilons = [np.inf]          # min distance in each generation
        self.generation_size = generation_size   # number of particles to be simulated at each generation
        self.all_simulated_data = []  # store all simulated data
        self.all_particles: List[candidateType] = []       # store all simulated particles
        self.all_distances: List[float] = []       # store all simulated distances
        self.maxiter = maxiter
        self.iterations = 0
        
    
    
    def calculate_distances_parallel(self,particles):
        with multiprocessing.Pool(self.cores) as p:
            res_map = p.map(self.simulator, particles)
        simulated_data = list(res_map)
        distances = [self.distance_function(self.Yobs, res) for res in simulated_data]

        # save all simulated results
        self.all_simulated_data.extend(simulated_data)
        self.all_distances.extend(distances)
        self.all_particles.extend(particles)
        
        return distances,simulated_data

    
    def check_t0_particles(self,particles):
        # check Topt < Tm, if false, resample from prior
        for particle in particles:
            for idp in particle.keys():
                if 'Tm' not in idp: continue
                id_topt = idp.split('_')[0]+'_Topt'
                if particle[idp]>particle[id_topt]: continue
                count = 0 # maximal resample times
                while count<10:
                    tm = self.posterior[idp].rvfv()
                    topt = self.posterior[id_topt].rvfv()
                    if tm>topt:
                        particle[idp] = tm
                        particle[id_topt]= topt
                        break
                    count += 1
        return particles
    
    def simulate_a_generation(self):
        particles_t, simulated_data_t, distances_t = [], [], []
        self.simulations += self.generation_size
        particles = [{idp: rv.rvfv() for idp,rv in self.posterior.items()} for _ in range(self.generation_size)]
        particles = self.check_t0_particles(particles)
        distances,simulated_data = self.calculate_distances_parallel(particles)
        
        particles_t.extend(particles)
        simulated_data_t.extend(simulated_data)
        distances_t.extend(distances)
        
        return particles_t, simulated_data_t, distances_t
    
    def update_population(self,particles_t, simulated_data_t, distances_t):
        logging.info('Updating population')
        # save first generation
        if len(self.population) == 0:
            self.population_t0 = particles_t
            self.distances_t0 = distances_t
            self.simulated_data_t0 = simulated_data_t
        
        
        combined_particles = np.array(self.population + particles_t)
        combined_distances = np.array(self.distances + distances_t)
        # We must be careful here as the results can be None due to hung-up processes
        combined_simulated = np.array(self.simulated_data + simulated_data_t, dtype=object)
        
        sort_index = np.argsort(combined_distances)
        self.population = list(combined_particles[sort_index][:self.population_size])
        self.distances = list(combined_distances[sort_index][:self.population_size])
        self.simulated_data = list(combined_simulated[sort_index][:self.population_size])
        self.epsilons.append(np.max(self.distances))
        
        logging.info(f"Model epsilon: {str(self.epsilons[-1])}")
        
        
    def update_posterior(self):
        logging.info('Updating posterior')
        parameters = dict()   # {'Protein_Tm':[]}
        for particle in self.population:
            for p,val in particle.items(): 
                lst = parameters.get(p,[])
                lst.append(val)
                parameters[p] =lst
        
        for p, lst in parameters.items():
            self.posterior[p] = RV('normal', loc = np.mean(lst), scale = np.std(lst))
        
    
    def run_simulation(self):
        self.iterations = len(self.epsilons)-1
        while self.iterations < self.maxiter:
            logging.info(f"Running iteration {self.iterations+1} of {self.maxiter}")
            if self.epsilons[-1] <= self.min_epsilon:
                logging.info("Bayesian fitting procedure ended successfully")
                break
            particles_t, simulated_data_t, distances_t = self.simulate_a_generation()
            self.update_population(particles_t, simulated_data_t, distances_t)
            self.update_posterior()
            pickle.dump(self,open(self.outfile,'wb'))
            self.iterations += 1
            #logging.info(f"epsilon: {self.epsilons[-1]}")
        else:
            logging.warning("Maximum number of iterations reached")
