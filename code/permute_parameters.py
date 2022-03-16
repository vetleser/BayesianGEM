import pandas as pd
import numpy as np

def permute_thermal_params(thermalparams: pd.DataFrame):
    permuted_params= thermalparams.copy()
    for param in ('dHTH', 'dSTS', 'dCpu'):
        permuted_params[param] = thermalparams[param].sample(frac=1)
    permuted_params[['Topt','dCpt']] = thermalparams[['Topt','dCpt']]
    return permuted_params


def permute_params(params: pd.DataFrame, rng: np.random.Generator):
    n_enzymes = params.shape[0]
    permuted_params = params.copy()
    index_order = rng.choice(n_enzymes, n_enzymes)
    # Shuffles the enzyme labels
    permuted_params.set_index(params.index[index_order])
    return permuted_params

    

