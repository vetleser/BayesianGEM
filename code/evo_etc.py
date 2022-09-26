#!/usr/bin/env python
# coding: utf-8

# #### Add constraint: Topt<Tm

# In[1]:


from abc import abstractmethod, abstractmethod
from itertools import repeat
import logging
from typing import Callable, Dict, Iterable, List, Set
import dill
import numpy as np
import numpy.typing as npt
import scipy.stats
import time
from multiprocessing import cpu_count
from pebble.concurrent.process import TimeoutError
import pebble
from random_sampler import RV
from abc import ABC

simResultType = Dict[str, npt.NDArray[np.float64]]
priorType = Dict[str, RV]
candidateType = Dict[str, float]
distanceArgType = Dict[str, npt.NDArray[np.float64]]



# In[2]:




class GA(ABC):
    def __init__(self,simulator: Callable[[candidateType], simResultType], priors : priorType, min_epsilon: float,
    distance_function: Callable[[distanceArgType, distanceArgType], float],
                 Yobs: distanceArgType, outfile: str,cores: int = cpu_count(),generation_size: int = 128, mutation_frequency: float = 1.,
                  selection_proportion: float = 0.5, maxiter: int = 100000, rng: np.random.Generator = None, n_children : int = 2, mutation_prob: float = 1.,
                  save_intermediate: bool = False):
        '''
        simulator:       a function that takes a dictionary of parameters as input. Ouput {'data':Ysim}
        priors:          a dictionary which use id of parameters as keys and RV class object as values
        min_epsilon:     minimal epsilon
        generation_size: the size of each population
        distance_function: a function that calculate the distance between observed data and simulated data
        Yobs:            observed data
        outfile:         unique id for the experiment. This will be also used to continue a simulation that 
                         is partly done
        cores:           number of treads
        mutation_frequency: The poisson lambda parameter to specify the number of mutations for each differentiation.
         Note that no mutation takes place in case the candidate is not selected for mutation (see parameter mutation_prob)
        selection_proportion: The proportion of the population in each generation used to produce offspring
        maxiter: The maximum number of generations to simulate
        n_children: The number of children to produce for each crossing
        mutation_prob: The probability that the mutator will apply mutation on a specific candidate solution
        save_intermediate: Should intermediate results be saved for each iteration? If toggled on, computations can be resumed if interrupted prematurly, but this will come at a performance penalty which parallelism cannot alleviate.

        
        !!!Important: distance is to be minimized!!!
        '''
        self.simulator = simulator
        self.priors = priors
        self.param_std: Dict[str, float] = {param: prior.scale for param, prior in priors.items()}
        if rng is None:
            default_seed = 1952
            self.rng = np.random.Generator(np.random.PCG64(default_seed))
        else:
            self.rng = rng
        self.distance_function = distance_function
        self.min_epsilon = min_epsilon
        self.Yobs = Yobs
        self.outfile = outfile
        # Compared to SMC-ABC this seems a bit odd. The reationale is indirection.
        # The indicies in the list specify which of the particles in self.all_particles is part of the current population
        self.population: List[List[int]] = []
        self.cores = cores    
        self.epsilons: List[float] = []          # min distance in each generation
        self.generation_size = generation_size   # number of particles to be simulated at each generation
        self.mutation_frequency = mutation_frequency
        self.selection_proportion = selection_proportion
        self.maxiter = maxiter
        # NOTE: all_simulated_data, all_distances, all_particles and birth_generation MUST be aligned
        self.all_simulated_data: List[simResultType] = []
        self.all_distances: List[float] = []
        self.all_particles: List[candidateType] = []
        # Only applicable for tournament replacers: The number of time particles has been in tournament
        # For truncation-like replacers: The number of times the particles has been selected at random without being one of the elites
        self.times_challenged: List[int] = []
        # Specifies the generation each of the particles are born in
        self.birth_generation: List[int] = []
        self.generation = 0
        self.n_children = n_children
        self.mutation_prob = mutation_prob
        self.save_intermediate = save_intermediate

        
    def evaluate_candiates(self, candidates: List[candidateType]):
        # Specifying timeout of 30 minutes
        timeout = 30*60
        # This function both evaluates newly born individuals and store them into the archive
        start = time.time()
        simulated_data = []
        try:
            # This code takes care of stalled parallel processes
            with pebble.ProcessPool(self.cores) as p:
                res_iter = p.map(self.simulator, candidates,timeout=timeout).result()
                while True:
                    try:
                        raw_res = res_iter.next()
                    except StopIteration:
                        # We have now iterated over all particles
                        break
                    except TimeoutError:
                        logging.info("Evaluation of particle time out")
                    else:
                        logging.info("Evaluation of particle ran successfully")
                        simulated_data.append(raw_res)
        except (OSError, RuntimeError) as e:
            logging.error('failed parallel_evaluation_mp: {0}'.format(str(e)))
            raise
        
        
        distances = [self.distance_function(self.Yobs, res) for res in simulated_data]

        # save all simulated results
        self.all_simulated_data.extend(simulated_data)
        self.all_distances.extend(distances)
        self.all_particles.extend(candidates)
        self.birth_generation.extend(repeat(self.generation,len(simulated_data)))
        self.times_challenged.extend(repeat(0,len(simulated_data)))
        end = time.time()
        logging.debug('Completed parallel evaluation of candiates in {0} seconds'.format(end - start))
        return


    def generator(self) -> candidateType:
            candidate = {param: prior.loc for param, prior in self.priors.items()}
            self.mutate(candidate=candidate)
            return candidate


    def mutate(self, candidate: candidateType):

        def mutate_param(entry: str):
            candidate[entry] = self.rng.normal(loc=candidate[entry], scale=self.param_std[entry])

        def correct_validity(entry: str) -> None:
            # As we only change one parameter at a time, we only need to check
            # the validity of the parameters of one enzyme
            protein_id = entry.split('_')[0]
            Tm_key = protein_id + "_Tm"
            Topt_key = protein_id + "_Topt"
            for _ in range(10):
                # We try to get things right 10 times before we give up
                Tm = candidate[Tm_key]
                Topt = candidate[Topt_key]
                if Tm < Topt:
                    mutate_param(Tm_key)
                    mutate_param(Topt_key)
                else:
                    break
        if float(self.rng.uniform(low=0, high=1)) < self.mutation_prob:
            n_mutations = self.rng.poisson(lam=self.mutation_frequency)
            entries: List[str] = list(candidate.keys())
            for _ in range(n_mutations):
                entry: str = self.rng.choice(entries)
                mutate_param(entry)
                correct_validity(entry)


    def create_offspring(self, mom: candidateType, dad: candidateType):
        children = self.cross(mom, dad)
        for child in children:
            self.mutate(child)
        return children



    def cross(self, mom: candidateType, dad: candidateType) -> Iterable[candidateType]:
        common_keys = set(mom.keys())
        common_keys.update(dad.keys())
        return [{key: (mom[key] + dad[key]) / 2 for key in common_keys} for _ in range(self.n_children)]

    def generate_children(self,parents: npt.NDArray[np.int64]):
        if len(parents) % 2 != 0:
                # If we have an add number of parents, remove the last one
                parents = parents[:-1]
        moms = parents[0::2]
        dads = parents[1::2]
        # Generate children
        logging.info(f"Generating {self.n_children*len(moms)} children")
        children = []
        for mom, dad in zip(moms, dads):
            children.extend(self.create_offspring(self.all_particles[mom],self.all_particles[dad]))
        logging.info(f"Evaluating fitness of children")
        self.evaluate_candiates(children)


    def update_std(self):
        """
        This method updates the enzyme parameter preceived standard deviations. Corresponds to update_posterior() in abc_etc 
        """
        logging.info('Updating standard deviations to parameters')
        parameters = dict()   # {'Protein_Tm':[]}
        for particle_idx in self.population[-1]:
            particle = self.all_particles[particle_idx]
            for p,val in particle.items(): 
                lst = parameters.get(p,[])
                lst.append(val)
                parameters[p] = lst
        
        for p, lst in parameters.items():
            self.param_std[p] = np.std(lst)

    
    def select_parents(self, current_population: npt.NDArray[np.int64]) -> npt.NDArray[np.int64]:
        population_fitness = np.array(self.all_distances)[current_population]
        n_parents = int(self.selection_proportion * self.generation_size)
        logging.info(f"Selecting {n_parents} parents to generating children")
        population_ranks = scipy.stats.rankdata(-population_fitness,method = "average") # Remember, this is a minimization problem, so lower values should have higher ranks
        parents = self.rng.choice(current_population, size=n_parents, replace=False, p=population_ranks / np.sum(population_ranks))
        return parents


    @abstractmethod
    def replacer(self, combined_population: npt.NDArray[np.int64]):
        pass


    def replace_population(self,combined_population: npt.NDArray[np.int64]):
        individuals_to_replace = len(combined_population) - self.generation_size
        logging.info(f"Replacing {individuals_to_replace} individuals in the population")
        # Remove individuals until population size is reached
        surviving_individuals = self.replacer(combined_population=combined_population)
        logging.info(f"Updating population")
        self.population.append(list(surviving_individuals))
            
    
    def particle_distance(self, idx_1: int, idx_2: int):
        epsilon = 1e-8
        particle_1 = self.all_particles[idx_1]
        particle_2 = self.all_particles[idx_2]
        return sum([(particle_1[key] - particle_2[key])**2 / (self.param_std[key]**2 + epsilon) for key in self.param_std.keys()])

    
    def simulate_generation(self):
        current_population = np.array(list(self.population[-1]))
        parents = self.select_parents(current_population=current_population)
        self.generate_children(parents=parents)
        children_idxs = np.flatnonzero(np.array(self.birth_generation) == self.generation)
        combined_population = np.concatenate([current_population,children_idxs])
        self.replace_population(combined_population)
        self.update_std()
        max_generation_epsilon = max(self.all_distances[p] for p in self.population[-1])
        self.epsilons.append(max_generation_epsilon)
        logging.info(f"Model epsilon {max_generation_epsilon}")
        

    
    def run_simulation(self) -> None:
        
        # Ensures random state is respected
        for item in self.priors.values():
            item.set_rng(rng=self.rng)
        
        if self.generation == 0:
            # Generating initial population
            logging.info("Generating initial population")
            initial_population = [self.generator() for _ in range(self.generation_size)]
            logging.info("Evaluating initial population")
            self.evaluate_candiates(initial_population)
            # Assign all individuals to be part of the first generation,
            # but beware, some of the evaluation tasks may have timed out
            self.population.append(list(range(len(self.all_particles))))
            self.update_std()
            max_generation_epsilon = max(self.all_distances[p] for p in self.population[-1])
            self.epsilons.append(max_generation_epsilon)
            logging.info(f"Model epsilon {max_generation_epsilon}")
            self.generation += 1
            if self.save_intermediate:
                    dill.dump(self,open(self.outfile,'wb'))
        
        
        while self.generation <= self.maxiter:
            if max_generation_epsilon < self.min_epsilon:
                logging.info(f"Fitness objective reached at generation {self.generation}")
                logging.info(f"Exiting evolution")
                break
            
            logging.info(f"Running generation {self.generation} of {self.maxiter}")
            self.simulate_generation()
            self.generation += 1
            if self.save_intermediate:
                    dill.dump(self,open(self.outfile,'wb'))
        else:
            # This else-clause belongs to the main evolution loop
            logging.info("Fitness objective not reached after maximum number of generations")
            logging.info("Exiting evolution")
        
        logging.info(f"Saving results to {self.outfile}")
        dill.dump(self, file=open(self.outfile,mode='wb'))



class TruncationGA(GA):
    def __init__(self, simulator: Callable[[candidateType], simResultType], priors: priorType, min_epsilon: float, distance_function: Callable[[distanceArgType, distanceArgType], float], Yobs: distanceArgType, outfile: str, cores: int = cpu_count(), generation_size: int = 128, mutation_frequency: float = 1, selection_proportion: float = 0.5, maxiter: int = 100000, rng: np.random.Generator = None, n_children: int = 2, mutation_prob: float = 1, save_intermediate: bool = False, num_elites: int=0):
        super().__init__(simulator, priors, min_epsilon, distance_function, Yobs, outfile, cores, generation_size, mutation_frequency, selection_proportion, maxiter, rng, n_children, mutation_prob, save_intermediate)
        self.num_elites = num_elites

    def replacer(self, combined_population: npt.NDArray[np.int64]):
        combined_population_fitness: npt.NDArray[np.int64] = np.array(self.all_distances)[combined_population]
        elites = combined_population[np.argsort(combined_population_fitness)[:self.num_elites]]
        non_elites = np.setdiff1d(combined_population,elites)
        lucky_freeloaders = self.rng.choice(non_elites,size= self.generation_size - len(elites))
        return np.concatenate((elites,lucky_freeloaders))




class TournamentGA(GA):
    def __init__(self, simulator: Callable[[candidateType], simResultType], priors: priorType, min_epsilon: float, distance_function: Callable[[distanceArgType, distanceArgType], float], Yobs: distanceArgType, outfile: str, cores: int = cpu_count(), generation_size: int = 128, mutation_frequency: float = 1, selection_proportion: float = 0.5, maxiter: int = 100000, rng: np.random.Generator = None, n_children: int = 2, mutation_prob: float = 1, save_intermediate: bool = False, locality: int=1):
        super().__init__(simulator, priors, min_epsilon, distance_function, Yobs, outfile, cores, generation_size, mutation_frequency, selection_proportion, maxiter, rng, n_children, mutation_prob, save_intermediate)
        self.locality = locality
    
    def replacer(self, combined_population: npt.NDArray[np.int64]):
        alive_individuals: Set[int] = set(combined_population)
        while len(alive_individuals) > self.generation_size:
            # Does a local tournament
            first = self.rng.choice(combined_population)
            if first not in alive_individuals:
                continue
            other_particles = np.array([other for other in alive_individuals if first != other])
            distance_to_first = np.array([self.particle_distance(first,other) for other in other_particles])
            particles_to_select = other_particles[np.argsort(distance_to_first)[:self.locality]]
            selected_particle = self.rng.choice(particles_to_select)

            if self.all_distances[first] > self.all_distances[selected_particle]:
                alive_individuals.remove(first)
            else:
                alive_individuals.remove(selected_particle)
        return alive_individuals

