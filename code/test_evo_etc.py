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
import GEMS
import os
import pandas as pd
import sys

from concurrent.futures import ProcessPoolExecutor, TimeoutError
import time

import dill

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s') 


logging.info('Imports done')


simResultType = Dict[str, npt.NDArray[np.float64]]
priorType = Dict[str, RV]
candidateType = Dict[str, float]
distanceArgType = Dict[str, npt.NDArray[np.float64]]

logging.info('Step 1')


generation_size = 128

path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)

logging.info('Step 2')

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

priors = dict()
for ind in params.index: 
    for col in ['Tm','Topt','dCpt']: 
        priors['{0}_{1}'.format(ind,col)] = RV('normal',
                                                    loc=params.loc[ind,col],
                                                    scale=params.loc[ind,col+'_std'])

parameter_names = list(priors.keys())

def generator() -> candidateType:
            candidate = {param: priors[param].rvfv() for param in parameter_names}
            #[self.correct_validity(candidate=candidate,entry=entry) for entry in candidate]
            return candidate

logging.info('Step 3')



cores: int = cpu_count()
simulator = GEMS.simulate_at_two_conditions_2
candidates = [generator() for _ in range(generation_size)]
timeout = 30*60

logging.info('Step 4')

logging.info('Trying to pickle candidates:')

try:
    dill.dumps(candidates)
    logging.info("Candidate is picklable.")
except Exception as e:
    logging.error(f"Candidate is not picklable: {e}")

logging.info('Trying to pickle simulator:')

try:
    dill.dumps(simulator)
    logging.info("simulator is picklable.")
except Exception as e:
    logging.error(f"Simulator is not picklable: {e}")


logging.info('Trying to wrap simulator:')

def wrapped_simulator(candidate):
    # Wrap the simulator with serialization/deserialization
    serialized_simulator = dill.dumps(simulator)
    simulator_func = dill.loads(serialized_simulator)
    return simulator_func(candidate)


for candidate in candidates:
    try:
        serialized = dill.dumps(wrapped_simulator(candidate))
        logging.info("Wrapped simulator works for candidate.")
    except Exception as e:
        logging.error(f"Serialization failed for candidate: {e}")
        break

# with ProcessPoolExecutor() as executor:
#     # Serialize candidates
#     serialized_candidates = [dill.dumps(candidate) for candidate in candidates]

#     # Use `executor.map` with the wrapped simulator
#     results = executor.map(wrapped_simulator, serialized_candidates)

#     # Deserialize results
#     res_iter = [dill.loads(res) for res in results]


try:
    # No need for creating a parallel cluster in this case
    # res_iter: Iterable[simResultType] = map(simulator,candidates)
    # while True:
    #     try:
    #         raw_res = res_iter.__next__()
    #     except StopIteration:
    #         # We have now iterated over all particles
    #         break
        
    with pebble.ProcessPool(cores) as p:
        res_iter: Iterable[simResultType] = p.map(wrapped_simulator, candidates, timeout=timeout).result()
        while True:
            try:
                raw_res = next(res_iter)
                # Process the result
            except StopIteration:
                logging.info("No more results available.")
                break
            except TimeoutError:
                logging.error("Timeout occurred while processing results.")
                continue  # Skip to the next result
            except Exception as e:
                logging.error(f"Unexpected error: {e}")
                p.stop()  # Explicitly stop the ProcessPool
                p.join()  # Wait for the pool to clean up
                sys.exit(1)  # Exit with failure code
except Exception as e:
    logging.error(f"Critical error in ProcessPool setup or execution: {e}")
    sys.exit(1)  # Exit with failure code

# with ProcessPoolExecutor(max_workers=cores) as executor:
#     futures = [executor.submit(simulator, candidate) for candidate in candidates]
#     try:
#         for future in futures:
#             raw_res = future.result(timeout=timeout)
#     except TimeoutError:
#         logging.error("Timeout occurred while processing results.")
#     except Exception as e:
#         logging.error(f"Unexpected error: {e}")

logging.info('Finished')

