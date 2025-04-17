#!/usr/bin/env python
# coding: utf-8

# #### Add constraint: Topt<Tm

# In[1]:

from itertools import repeat
import logging
from typing import Callable, Dict, Iterable, List, Set
import dill
import numpy as np
import numpy.typing as npt
import time
from multiprocessing import cpu_count
from pebble.concurrent.process import TimeoutError
import pebble
from random_sampler import RV
import copy

simResultType = Dict[str, npt.NDArray[np.float64]]
priorType = Dict[str, RV]
candidateType = Dict[str, float]
distanceArgType = Dict[str, npt.NDArray[np.float64]]


class SimulatedAnnealing():
    def __init__(self, 
                 simulator: Callable[[candidateType], simResultType], 
                 priors : priorType, 
                 min_epsilon: float,
                 distance_function: Callable[[distanceArgType, distanceArgType], float],
                 Yobs: distanceArgType, 
                 outfile: str,
                 cores: int = cpu_count(),
                 maxiter: int = 100000, 
                 rng: np.random.Generator = None, 
                 save_intermediate: bool = False,
                 
                 generation_size: int = 128,
                 cooling_rate: float = 0.95,
                 initial_temp: float = 100,
                 final_temp: float = 0.1 
                 ):
        """Implements the Simulated Annealing algorithm designed to detect multiple optima in the fitness landscape

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
            
            !!! Distance is to be minimized!!!
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
        self.save_intermediate = save_intermediate

        self.cooling_rate = cooling_rate
        self.initial_temp = initial_temp
        self.current_temp = initial_temp
        self.final_temp = final_temp


    def generator(self) -> candidateType:
            candidate = {param: self.priors[param].rvfv() for param in self.parameter_names}
            [self.correct_validity(candidate=candidate,entry=entry) for entry in candidate]
            return candidate


    def correct_validity(self, candidate, entry: str) -> None:
        # As we only change one parameter at a time, we only need to check
        # the validity of the parameters of one enzyme
        # Topt > Tm > 0 (in Kelvin of course) in real life, but mutation may disregard this constraint, so we have to account for it
        split_entry = entry.split('_')
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
            if not Tm > Topt > 0:
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


    def mutate_param(self,candidate, entry: str):
            candidate[entry] = self.priors[entry].rvfv()

    def check_validity(self, candidate, entry: str) -> bool:
        # As we only change one parameter at a time, we only need to check
        # the validity of the parameters of one enzyme
        # Topt > Tm in real life, but mutation may disregard this constraint, so we have to account for it
        split_entry = entry.split('_')
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
        return Tm > Topt > 0


    def particle_distance(self, idx_1: int, idx_2: int)-> float:
        epsilon = 1e-8
        particle_1 = self.all_particles[idx_1]
        particle_2 = self.all_particles[idx_2]
        return sum([(particle_1[key] - particle_2[key])**2 / (self.param_std[key]**2 + epsilon) for key in self.param_std.keys()])



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


    def change_all_parameters(self, old_particle):
        logging.info("Changing all parameters of a particle")
        new_particle = copy.deepcopy(self.all_particles[old_particle])
        for key in new_particle:
            old_parameter_value = new_particle[key]
            new_particle[key] += np.random.normal(0, 0.1)
            if not self.check_validity(new_particle,key):
                # This new parameter value violates our constraints, so we must revert the change
                new_particle[key] = old_parameter_value
        return new_particle

    def simple_evaluation(self, particle1, particle2):
        return particle2
        


    def update_population(self, original_population: npt.NDArray[np.int64]):
        logging.info(f"Applying Simulated Annealing to current population")
        new_population: Set[int] = set(original_population)
        for particle in original_population:
            new_particle = self.change_all_parameters(particle)
            chosen_particle = self.simple_evaluation(particle, new_particle)
            new_population.remove(particle)
            new_population.add(chosen_particle)
        self.population.append(list(new_population))
        return

        

   

    def simulate_generation(self):
        current_population = np.array(list(self.population[-1]))

        #Remove
        # primary_parents = self.select_parents(current_population=current_population)
        # secondary_parents = np.setdiff1d(current_population,primary_parents)
        # self.generate_children(primary_parents,secondary_parents)
        # children_idxs = np.flatnonzero(np.array(self.birth_generation) == self.generation)
        
        
        
        #Replace
        #self.replace_population(current_population,children=children_idxs)
        self.update_population(current_population)

        max_generation_epsilon = max(self.all_distances[p] for p in self.population[-1])
        self.epsilons.append(max_generation_epsilon)
        self.update_std()
        #logging.info(f"Model epsilon {max_generation_epsilon}")


    def run_simulation(self) -> None:
        # Ensures random state is respected
        for item in self.priors.values():
            item.set_rng(rng=self.rng)
        
        if self.generation == 0:
            # Generating initial population
            logging.info("Generating initial population")
            initial_population = [self.generator() for _ in range(self.generation_size)]
            logging.info("Evaluating initial population")
            #self.evaluate_candiates(initial_population)
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
            if self.current_temp < self.final_temp:
                logging.info(f"Temperature has reached a minimum at generation {self.generation}")
            
            logging.info(f"Running generation {self.generation} of {self.maxiter}. Current temperature: {self.current_temp}")
            self.simulate_generation()
            self.generation += 1
            self.current_temp *= self.cooling_rate
            if self.save_intermediate:
                    dill.dump(self,open(self.outfile,'wb'))
        else:
            # This else-clause belongs to the main evolution loop
            logging.info("Fitness objective not reached after maximum number of generations")
            logging.info("Exiting evolution")
        
        logging.info(f"Saving results to {self.outfile}")
        dill.dump(self, file=open(self.outfile,mode='wb'))
