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
        # For tournament replacers: The number of time particles has been in tournament
        # For truncation-like replacers: The number of times the particles has been alive and not being one of the elites
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
        # Candidates for which evaluating fitness was successful
        successfull_candiates = set()
        candidate_counter = 0
        if self.cores == 1:
            # No need for creating a parallel cluster in this case
            res_iter: Iterable[simResultType] = map(self.simulator,candidates)
            while True:
                try:
                    raw_res = res_iter.__next__()
                except StopIteration:
                    # We have now iterated over all particles
                    break
                else:
                    logging.info("Evaluation of particle ran successfully")
                    simulated_data.append(raw_res)
                    successfull_candiates.add(candidate_counter)
                finally:
                    candidate_counter += 1
        else:
                try:
                    # This code takes care of stalled parallel processes
                    with pebble.ProcessPool(self.cores) as p:
                        res_iter: Iterable[simResultType] = p.map(self.simulator, candidates,timeout=timeout).result()
                        while True:
                            try:
                                raw_res = res_iter.next()
                            except StopIteration:
                                # We have now iterated over all particles
                                break
                            except TimeoutError:
                                logging.info("Evaluation of particle timed out")
                            else:
                                logging.info("Evaluation of particle ran successfully")
                                simulated_data.append(raw_res)
                                successfull_candiates.add(candidate_counter)
                            finally:
                                candidate_counter += 1
                except (OSError, RuntimeError) as e:
                    logging.error('failed parallel_evaluation_mp: {0}'.format(str(e)))
                    raise
        
        
        distances = [self.distance_function(self.Yobs, res) for res in simulated_data]

        # save all simulated results
        self.all_simulated_data.extend(simulated_data)
        self.all_distances.extend(distances)
        # This deals with the problem of candidates failing evaluation
        self.all_particles.extend([candidate for counter, candidate in enumerate(candidates) if counter in successfull_candiates])
        self.birth_generation.extend(repeat(self.generation,len(simulated_data)))
        self.times_challenged.extend(repeat(0,len(simulated_data)))
        end = time.time()
        logging.debug('Completed parallel evaluation of candiates in {0} seconds'.format(end - start))
        return


    def generator(self) -> candidateType:
            candidate = {param: prior.rvfv() for param, prior in self.priors.items()}
            [self.correct_validity(candidate=candidate,entry=entry) for entry in candidate]
            return candidate

    def correct_validity(self, candidate, entry: str) -> None:
        # As we only change one parameter at a time, we only need to check
        # the validity of the parameters of one enzyme
        # Topt > Tm in real life, but mutation may disregard this constraint, so we have to account for it
        split_entry = entry.split('_')[0]
        # We assume that entries are of the form PROTID_{Tm,Topt,dCpt}
        # If this is not the case, we assume that the algorithm is used for another kind of inference problem,
        # so we skip this domain-specific check. This also applies to the dCPt as mutatating them does not violate the constraint
        if len(split_entry) != 2 or split_entry[1] not in ("Tm","Topt"):
            return
        protein_id = split_entry[0]
        Tm_key = protein_id + "_Tm"
        Topt_key = protein_id + "_Topt"
        for _ in range(10):
            # We try to get things right 10 times before we give up
            Tm = candidate[Tm_key]
            Topt = candidate[Topt_key]
            if Tm < Topt:
                self.mutate_param(candidate,Tm_key)
                self.mutate_param(candidate,Topt_key)
            else:
                break
        
    def mutate_param(self,candidate, entry: str):
            candidate[entry] = self.rng.normal(loc=candidate[entry], scale=self.param_std[entry])

    def mutate(self, candidate: candidateType):
        if float(self.rng.uniform(low=0, high=1)) < self.mutation_prob:
            n_mutations = self.rng.poisson(lam=self.mutation_frequency)
            entries: List[str] = list(candidate.keys())
            for _ in range(n_mutations):
                entry: str = self.rng.choice(entries)
                self.mutate_param(candidate,entry)
                self.correct_validity(candidate,entry)


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
    def replacer(self, combined_population: npt.NDArray[np.int64])-> List[int]:
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
        for non_elite in non_elites:
            self.times_challenged[non_elite] +=1
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
            self.times_challenged[first] += 1
            self.times_challenged[selected_particle] += 1

            if self.all_distances[first] > self.all_distances[selected_particle]:
                alive_individuals.remove(first)
            else:
                alive_individuals.remove(selected_particle)
        return alive_individuals


class UniformParentTournamentGA(TournamentGA):
    def __init__(self, simulator: Callable[[candidateType], simResultType], priors: priorType, min_epsilon: float, distance_function: Callable[[distanceArgType, distanceArgType], float], Yobs: distanceArgType, outfile: str, cores: int = cpu_count(), generation_size: int = 128, mutation_frequency: float = 1, selection_proportion: float = 0.5, maxiter: int = 100000, rng: np.random.Generator = None, n_children: int = 2, mutation_prob: float = 1, save_intermediate: bool = False, locality: int = 1):
        super().__init__(simulator, priors, min_epsilon, distance_function, Yobs, outfile, cores, generation_size, mutation_frequency, selection_proportion, maxiter, rng, n_children, mutation_prob, save_intermediate, locality)
    def select_parents(self, current_population: npt.NDArray[np.int64]) -> npt.NDArray[np.int64]:
        n_parents = int(self.selection_proportion * self.generation_size)
        logging.info(f"Selecting {n_parents} parents to generating children")
        parents = self.rng.choice(current_population, size=n_parents, replace=False)
        return parents


class CrowdingDE():
    def __init__(self, simulator: Callable[[candidateType], simResultType], priors : priorType, min_epsilon: float,
    distance_function: Callable[[distanceArgType, distanceArgType], float],
                 Yobs: distanceArgType, outfile: str,cores: int = cpu_count(),
                 generation_size: int = 128, scaling_factor: float = 0.5, crossover_prob = 0.9,
                 maxiter: int = 100000, rng: np.random.Generator = None, n_children : int = 2,
                  save_intermediate: bool = False):
        """Implements the Crowding Differential Evolution algorithm designed to detect multiple optima in the fitness landscape

        Args:
            simulator (Callable[[candidateType], simResultType]): a function that takes a dictionary of parameters as input. Ouput {'data':Ysim}
            priors (priorType): a dictionary which use id of parameters as keys and RV class object as values
            min_epsilon (float): minimal epsilon
            distance_function (Callable[[distanceArgType, distanceArgType], float]): a function that calculate the distance between observed data and simulated data
            Yobs (distanceArgType): observed data
            outfile (str): unique id for the experiment. This will be also used to continue a simulation that 
                         is partly done
            cores (int, optional): number of treads. Defaults to cpu_count().
            generation_size (int, optional): the size of each population. Defaults to 128.
            maxiter (int, optional): The maximum number of generations to simulate. Defaults to 100000.
            rng (np.random.Generator, optional): Random number generator in order to ensure reproducibility. Defaults to None.
            n_children (int, optional): The number of children of generate for each generation. Defaults to 2.
            scaling_factor (float, optional): The scaling factor for the Differential Evolution, meaning how much weight is applied to the two secondary individuals. Defaults to 0.5.
            crossover_prob (float, optional): The crossover probability for the two secondary parents. Defaults to 0.9.
            save_intermediate (bool, optional): Should intermediate results be saved for each iteration? If toggled on, computations can be resumed if interrupted prematurly, but this will come at a performance penalty which parallelism cannot alleviate.. Defaults to False.
        """
        self.simulator = simulator
        self.priors = priors
        self.parameter_names = list(priors.keys())
        self.n_parameters = len(self.parameter_names)
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
        self.maxiter = maxiter
        # NOTE: all_simulated_data, all_distances, all_particles and birth_generation MUST be aligned
        self.all_simulated_data: List[simResultType] = []
        self.all_distances: List[float] = []
        self.all_particles: List[candidateType] = []
        # The number of times the particle has been challenged by an offspring
        self.times_challenged: List[int] = []
        # Specifies the generation each of the particles are born in
        self.birth_generation: List[int] = []
        self.generation = 0
        self.n_children = n_children
        self.crossover_prob = crossover_prob
        self.scaling_factor = scaling_factor
        self.save_intermediate = save_intermediate


    def generator(self) -> candidateType:
            candidate = {param: prior.rvfv() for param, prior in self.priors.items()}
            [self.correct_validity(candidate=candidate,entry=entry) for entry in candidate]
            return candidate


    def correct_validity(self, candidate, entry: str) -> None:
        # As we only change one parameter at a time, we only need to check
        # the validity of the parameters of one enzyme
        # Topt > Tm > 0 (in Kelvin of course) in real life, but mutation may disregard this constraint, so we have to account for it
        split_entry = entry.split('_')[0]
        # We assume that entries are of the form PROTID_{Tm,Topt,dCpt}
        # If this is not the case, we assume that the algorithm is used for another kind of inference problem,
        # so we skip this domain-specific check. This also applies to the dCPt as mutatating them does not violate the constraint
        if len(split_entry) != 2 or split_entry[1] not in ("Tm","Topt"):
            return
        protein_id = split_entry[0]
        Tm_key = protein_id + "_Tm"
        Topt_key = protein_id + "_Topt"
        for _ in range(10):
            # We try to get things right 10 times before we give up
            Tm = candidate[Tm_key]
            Topt = candidate[Topt_key]
            if Tm < Topt or Tm < 0:
                self.mutate_param(candidate,Tm_key)
                self.mutate_param(candidate,Topt_key)
            else:
                break
    
    
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

    def evaluate_candiates(self, candidates: List[candidateType]):
        # Specifying timeout of 30 minutes
        timeout = 30*60
        # This function both evaluates newly born individuals and store them into the archive
        start = time.time()
        simulated_data = []
        # Candidates for which evaluating fitness was successful
        successfull_candiates = set()
        candidate_counter = 0
        if self.cores == 1:
            # No need for creating a parallel cluster in this case
            res_iter: Iterable[simResultType] = map(self.simulator,candidates)
            while True:
                try:
                    raw_res = res_iter.__next__()
                except StopIteration:
                    # We have now iterated over all particles
                    break
                else:
                    logging.info("Evaluation of particle ran successfully")
                    simulated_data.append(raw_res)
                    successfull_candiates.add(candidate_counter)
                finally:
                    candidate_counter += 1
        else:
                try:
                    # This code takes care of stalled parallel processes
                    with pebble.ProcessPool(self.cores) as p:
                        res_iter: Iterable[simResultType] = p.map(self.simulator, candidates,timeout=timeout).result()
                        while True:
                            try:
                                raw_res = res_iter.next()
                            except StopIteration:
                                # We have now iterated over all particles
                                break
                            except TimeoutError:
                                logging.info("Evaluation of particle timed out")
                            else:
                                logging.info("Evaluation of particle ran successfully")
                                simulated_data.append(raw_res)
                                successfull_candiates.add(candidate_counter)
                            finally:
                                candidate_counter += 1
                except (OSError, RuntimeError) as e:
                    logging.error('failed parallel_evaluation_mp: {0}'.format(str(e)))
                    raise
        
        distances = [self.distance_function(self.Yobs, res) for res in simulated_data]

        # save all simulated results
        self.all_simulated_data.extend(simulated_data)
        self.all_distances.extend(distances)
        # This deals with the problem of candidates failing evaluation
        self.all_particles.extend([candidate for counter, candidate in enumerate(candidates) if counter in successfull_candiates])
        self.birth_generation.extend(repeat(self.generation,len(simulated_data)))
        self.times_challenged.extend(repeat(0,len(simulated_data)))
        end = time.time()
        logging.debug('Completed parallel evaluation of candiates in {0} seconds'.format(end - start))
        return

    def select_parents(self, current_population: npt.NDArray[np.int64]) -> npt.NDArray[np.int64]:
            parents = self.rng.choice(current_population, size=self.n_children, replace=False)
            return parents

    def mutate_param(self,candidate, entry: str):
            candidate[entry] = self.priors[entry].rvf()

    def check_validity(self, candidate, entry: str) -> bool:
        # As we only change one parameter at a time, we only need to check
        # the validity of the parameters of one enzyme
        # Topt > Tm in real life, but mutation may disregard this constraint, so we have to account for it
        split_entry = entry.split('_')[0]
        # We assume that entries are of the form PROTID_{Tm,Topt,dCpt}
        # If this is not the case, we assume that the algorithm is used for another kind of inference problem,
        # so we skip this domain-specific check. This also applies to the dCPt as mutatating them does not violate the constraint
        if len(split_entry) != 2 or split_entry[1] not in ("Tm","Topt"):
            return True
        protein_id = split_entry[0]
        Tm_key = protein_id + "_Tm"
        Topt_key = protein_id + "_Topt"
        Tm = candidate[Tm_key]
        Topt = candidate[Topt_key]
        return Topt < Tm < 0


    def particle_distance(self, idx_1: int, idx_2: int)-> float:
        epsilon = 1e-8
        particle_1 = self.all_particles[idx_1]
        particle_2 = self.all_particles[idx_2]
        return sum([(particle_1[key] - particle_2[key])**2 / (self.param_std[key]**2 + epsilon) for key in self.param_std.keys()])


    def create_offspring(self,primary_parent: int,parent_2: int,parent_3: int)-> candidateType:
        parent_2_particle = self.all_particles[parent_2]
        parent_3_particle = self.all_particles[parent_3]
        offspring_particle = {parameter: value for parameter, value in self.all_particles[primary_parent].items()}
        mutate_start = self.rng.choice(self.n_parameters)
        for i in range(mutate_start,mutate_start+len(self.parameter_names)):
            if self.rng.uniform() > self.crossover_prob:
                break
            parameter_to_mutate = self.parameter_names[i % self.n_parameters]
            old_parameter_value = offspring_particle[parameter_to_mutate]
            offspring_particle[parameter_to_mutate] += self.scaling_factor*(parent_2_particle[parameter_to_mutate]-parent_3_particle[parameter_to_mutate])
            if not self.check_validity(offspring_particle,parameter_to_mutate):
                # This new parameter value violates our constraints, so we must revert the change
                offspring_particle[parameter_to_mutate] = old_parameter_value
        return offspring_particle

    def replace_population(self,original_population: npt.NDArray[np.int64], children: npt.NDArray[np.int64]):
        logging.info(f"Replacing population with children")
        # Remove individuals until population size is reached
        current_population: Set[int] = set(original_population)
        for child in children:
            other_particles = np.fromiter(current_population,int,len(current_population))
            a = [self.particle_distance(child,other) for other in other_particles]
            distance_to_child: npt.NDArray[np.float64] = np.array([self.particle_distance(child,other) for other in other_particles])
            closest_to_child = other_particles[np.argmin(distance_to_child)]
            if self.all_distances[closest_to_child] > self.all_distances[child]:
                # The child replaces an old individual
                current_population.remove(closest_to_child)
                current_population.add(child)
        logging.info(f"Updating population")
        self.population.append(list(current_population))

    def generate_children(self, primary_parents: npt.NDArray[np.int64],secondary_parents: npt.NDArray[np.int64]):
        # Generate children
        logging.info(f"Generating {self.n_children} children")
        children: List[candidateType] = []
        for primary_parent in primary_parents:
            parent_2, parent_3 = self.rng.choice(secondary_parents,size=2,replace=False)
            children.append(self.create_offspring(primary_parent, parent_2, parent_3))
        logging.info(f"Evaluating fitness of children")
        self.evaluate_candiates(children)

    def simulate_generation(self):
        current_population = np.array(list(self.population[-1]))
        primary_parents = self.select_parents(current_population=current_population)
        secondary_parents = np.setdiff1d(current_population,primary_parents)
        self.generate_children(primary_parents,secondary_parents)
        children_idxs = np.flatnonzero(np.array(self.birth_generation) == self.generation)
        self.replace_population(current_population,children=children_idxs)
        max_generation_epsilon = max(self.all_distances[p] for p in self.population[-1])
        self.epsilons.append(max_generation_epsilon)
        self.update_std()
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
            max_generation_epsilon = max(self.all_distances[p] for p in self.population[-1])
            self.epsilons.append(max_generation_epsilon)
            self.update_std()
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
