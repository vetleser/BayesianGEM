#!/usr/bin/env python
# coding: utf-8

# #### Add constraint: Topt<Tm

# In[1]:


import logging
from typing import Callable, Dict, Iterable, List
import dill
from inspyred.ec import replacers, variators
from inspyred.ec.variators import mutators
import numpy as np
import numpy.typing as npt
import time
import scipy.stats as ss
from scipy.stats import poisson
from multiprocessing import cpu_count
from random_sampler import RV
import inspyred
import inspyred.ec
import inspyred.ec.variators as variators
import inspyred.ec.selectors as selectors
import pathos

simResultType = Dict[str, Dict[str, npt.NDArray[np.float64]]]
priorType = Dict[str, RV]
candidateType = Dict[str, float]
distanceArgType = Dict[str, npt.NDArray[np.float64]]


# In[]

# This is a really ugly trick to use a numpy RNG
# when inspyred expects the generator to be of
# type random.Random
# For this to work, we must alias random.Random.sample
# using np.random.Generator.choice
class EvolutionGenerator(np.random.Generator):

    def __init__(self, bit_generator: np.random.BitGenerator) -> None:
        super().__init__(bit_generator)
    
    def sample(self, population, k, *, counts=None):
        return self.choice(a=population, size=k)

    @classmethod
    def from_numpy_generator(cls, generator: np.random.Generator):
        return cls(generator.bit_generator)



def default_rng(seed=None):
    return EvolutionGenerator(np.random.PCG64(seed))

# In[2]:




class GA:
    def __init__(self,simulator: Callable[[candidateType], simResultType], priors : priorType, min_epsilon: float,
    distance_function: Callable[[distanceArgType, distanceArgType], float],
                 Yobs: distanceArgType, outfile: str,cores: int = cpu_count(),generation_size: int = 128, mutation_frequency: float = 1.,
                  selection_proportion: float = 0.5, maxiter: int = 100000, rng: EvolutionGenerator = None, n_children : int = 2, mutation_prob: float = 1.,
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
            self.rng = default_rng(default_seed)
        else:
            self.rng = rng
        # Ensures random state is respected
        for item in self.priors.values():
            item.set_rng(rng=rng)
        self.distance_function = distance_function
        self.min_epsilon = min_epsilon
        self.Yobs = Yobs
        self.outfile = outfile
        self.population: List[inspyred.ec.Individual] = []
        self.cores = cores    
        self.epsilons = [np.inf]          # min distance in each generation
        self.generation_size = generation_size   # number of particles to be simulated at each generation
        self.mutation_frequency = mutation_frequency
        self.selection_proportion = selection_proportion
        self.maxiter = maxiter
        self.all_simulated_data = []
        self.all_distances = []
        self.all_particles = []
        self.generations = 0
        self.n_children = n_children
        self.mutation_prob = mutation_prob
        self.save_intermediate = save_intermediate

        

        def archiver(random, population, archive, args):
            # We do not really care about the inspyred archive, but
            # have to handle it anyway
            if archive is None:
                archive = []
            else:
                archive.append(None)
            self.population = population
            max_generation_epsilon = max(p.fitness for p in population)
            self.epsilons.append(max_generation_epsilon)
            logging.info(f"Model epsilon {max_generation_epsilon}")
            self.generations += 1

            if self.save_intermediate:
                    dill.dump(self,open(self.outfile,'wb'))
            return archive


        def terminator(population, num_generations, num_evaluations, args) -> bool:
            n_generations = self.generations
            logging.info(f"Running generation {n_generations+1} of {self.maxiter}")
            return self.epsilons[-1] <= self.min_epsilon or self.generations > self.maxiter

        # @evaluators.evaluator
        """ def evaluator(candidate: individualType, args: dict) -> float:
            res = self.simulator(args = {key: value.loc for key, value in candidate.items()})
            distance = self.distance_function(Yobs, res)
            return res, distance """

        def parallel_evaluation_mp(candidates: List[candidateType], args: dict):
            # Copy-catted from https://github.com/aarongarrett/inspyred/blob/master/inspyred/ec/evaluators.py
            logger = args['_ec'].logger
            
            try:
                nprocs: int = args['mp_nprocs']
            except KeyError:
                nprocs = pathos.multiprocessing.cpu_count()
                
            pickled_args = {}
            for key in args:
                try:
                    dill.dumps(args[key])
                    pickled_args[key] = args[key]
                except (TypeError, dill.PicklingError):
                    logger.debug('unable to pickle args parameter {0} in parallel_evaluation_mp'.format(key))
                    pass

            start = time.time()
            try:
                # Use parallel backend such as in abc_etc.py
                Q = pathos.helpers.mp.Manager().Queue()
                jobs = [pathos.helpers.mp.Process(target=self.simulate_one,args=(particle,index,Q)) 
                               for index,particle in enumerate(candidates)]
                
                for p in jobs: p.start()
                # The timeout is an emergency hatch designed to catch 
                # processes which for some reason are caught in a deadlock
                for p in jobs: p.join(timeout=1000)

                # Q may not always contain the result of all jobs we passed to it,
                # this must be handled carefully
                # We solve the problem by assuming that any jobs which fails to complete corresponds to zero fitness
                result_list = [None for _ in range(len(candidates))]

                while not Q.empty():
                    index,res = Q.get(timeout=1)
                    result_list[index] = res
                
            except (OSError, RuntimeError) as e:
                logger.error('failed parallel_evaluation_mp: {0}'.format(str(e)))
                raise
            else:
                all_simulated_data = []
                all_distances = []
                for entry in result_list:
                    if entry is None:
                        all_simulated_data.append(None)
                        all_distances.append(np.inf)
                    else:
                        res = entry 
                        distance = self.distance_function(self.Yobs,res)
                        all_simulated_data.append(res)
                        all_distances.append(distance)
                self.all_distances.extend(all_distances)
                self.all_simulated_data.extend(all_simulated_data)
                self.all_particles.extend(candidates)

                end = time.time()
                logger.debug('completed parallel_evaluation_mp in {0} seconds'.format(end - start))
                return all_distances


        def generator(random: EvolutionGenerator, args) -> candidateType:
            # def correct_validity(seed: individualType):
            #     seed = deepcopy(self.priors)
            #     # check Topt < Tm, if false, resample from prior
            #     protein_ids = set(entry.split('_')[0] for entry in seed)
            #     for id in protein_ids:
            #         Tm_key = id + "_Tm"
            #         Topt_key = id + "_Topt"
            #         for _ in range(10):
            #             # We try to get things right 10 times before we give up
            #             Tm = seed[Tm_key]
            #             Topt = seed[Topt_key]
            #             if Tm < Topt:
            #                 seed[Tm_key] = Tm.mutate()
            #                 seed[Topt_key] = Topt.mutate()
            #             else:
            #                 break
            #     return seed
                candidate = {param: prior.loc for param, prior in self.priors.items()}
                return mutate(random = random, candidate=candidate, args=args)

        # @variators.mutator

        def mutate(random: EvolutionGenerator, candidate: candidateType, args: dict) -> candidateType:

            def mutate_param(entry: str):
                candidate[entry] = random.normal(loc=candidate[entry], scale=self.param_std[entry])

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
            if float(random.uniform(low=0, high=1)) < self.mutation_prob:
                n_mutations = random.poisson(lam=self.mutation_frequency)
                entries: list[str] = list(candidate.keys())
                for _ in range(n_mutations):
                    entry: str = random.choice(entries)
                    mutate_param(entry)
                    correct_validity(entry)
            return candidate


        @variators.crossover
        def cross(random: EvolutionGenerator, mom: candidateType, dad: candidateType, args: dict) -> Iterable[candidateType]:
            common_keys = set(mom.keys())
            common_keys.update(dad.keys())
            return ({key: (mom[key] + dad[key]) / 2 for key in common_keys} for _ in range(self.n_children))


        self.archiver = archiver
        self.terminator = terminator
        # self.evaluator = evaluator
        self.generator = generator
        self.mutator = mutators.mutator(mutate)
        self.crossover = cross
        self.parallel_evaluation_mp = parallel_evaluation_mp
    

    def update_std(self):
        """
        This method updates the enzyme parameter preceived standard deviations. Corresponds to update_posterior() in abc_etc 
        """
        logging.info('Updating standard deviations to parameters')
        parameters = dict()   # {'Protein_Tm':[]}
        for particle in self.population:
            for p,val in particle.items(): 
                lst = parameters.get(p,[])
                lst.append(val)
                parameters[p] = lst
        
        for p, lst in parameters.items():
            self.param_std[p] = np.std(lst)



    def simulate_one(self,particle,index,Q):
        '''
        particle:  parameters 
        Q:      a multiprocessing.Queue object
        index:  the index in particles list
        '''

        res = self.simulator(particle)
        # ysim = {simulated}

        Q.put((index,res))

    
    def run_simulation(self) -> None:
        algorithm = inspyred.ec.EvolutionaryComputation(self.rng)
        algorithm.terminator = self.terminator
        algorithm.archiver = self.archiver
        algorithm.variator = [self.crossover, self.mutator]
        algorithm.selector = selectors.tournament_selection
        algorithm.replacer = replacers.truncation_replacement
        # Restore rng. This is a workaround due to https://github.com/uqfoundation/dill/issues/442
        self.rng = EvolutionGenerator.from_numpy_generator(self.rng)
        if self.generations > 0:
            # Restore computations from stored results
            seeds : List[candidateType] = [individual.candidate for individual  in self.population]
        else:
            seeds = None


        kwargsdict = {'mp_nprocs': self.cores, 'num_selected': int(self.generation_size * self.selection_proportion)}
        algorithm.evolve(generator=self.generator, evaluator=self.parallel_evaluation_mp, maximize=False, seeds = seeds, pop_size=self.generation_size, **kwargsdict)
        dill.dump(self, file=open(self.outfile,mode='wb'))
