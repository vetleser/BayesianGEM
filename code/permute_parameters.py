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
    Topt_order = rng.choice(n_enzymes, n_enzymes)
    permuted_params["Topt"] = params["Topt"].to_numpy()[Topt_order]
    permuted_params["Topt_std"] = params["Topt_std"].to_numpy()[Topt_order]
    Tm_order = rng.choice(n_enzymes, n_enzymes)
    permuted_params["Tm"] = params["Tm"].to_numpy()[Tm_order]
    permuted_params["Tm_std"] = params["Tm_std"].to_numpy()[Tm_order]
    dCpt_order = rng.choice(n_enzymes, n_enzymes)
    permuted_params["dCpt"] = params["dCpt"].to_numpy()[dCpt_order]
    permuted_params["dCpt_std"] = params["dCpt_std"].to_numpy()[dCpt_order]
    T90_order = rng.choice(n_enzymes, n_enzymes)
    permuted_params["T90"] = params["T90"].to_numpy()[T90_order]
    Length_order = rng.choice(n_enzymes, n_enzymes)
    permuted_params["Length"] = params["Length"].to_numpy()[Length_order]
    permuted_params.set_index(params.index)
    return permuted_params

    

