# In[]
from typing import Dict
import logging
import GEMS
import os
import pandas as pd
from etcpy import etc
import pickle
import GEMS
import numpy as np
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

candidateType = Dict[str, float]


def aerobic_fva(thermalParams: candidateType):
    """
    Run FVA on aerobic conditions

    Args:
        thermalParams: A dictionary of the model's thermal parameters
    """
    # thermalParams: a dictionary with ids like uniprotid_Topt
    param_dict = GEMS.format_input(thermalParams)
    mae = pickle.load(open(os.path.join(path,'models/aerobic.pkl'),'rb'))
    rae = etc.simulate_fva(mae,dfae_batch.index+273.15,param_dict=param_dict,sigma=0.5)
    rae["dataset"] = "aerobic"
    return rae


def anerobic_reduced_fva(thermalParams: candidateType):
    """
    Run FVA under anaerobic conditions

    Args:
        thermalParams: A dictionary of the model's thermal parameters
    """
    # thermalParams: a dictionary with ids like uniprotid_Topt 
    param_dict = GEMS.format_input(thermalParams)
    mae = pickle.load(open(os.path.join(path,'models/anaerobic.pkl'),'rb'))
    sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]
    ran = etc.simulate_fva(mae,np.array(sel_temp+273.15,param_dict=param_dict,sigma=0.5))
    ran["dataset"] = "anaerobic"
    return ran


def chemostat_fva(thermalParams):
    """
    Run FVA under chemostat conditions

    Args:
        thermalParams: A dictionary of the model's thermal parameters
    """
    param_dict = GEMS.format_input(thermalParams)
    mae = pickle.load(open(os.path.join(path,'models/aerobic.pkl'),'rb'))
    growth_id = 'r_2111'
    glc_up_id = 'r_1714_REV'
    prot_pool_id = 'prot_pool_exchange'
    dilut = 0.1
    sigma = 0.5
    
    solution = etc.fva_chemostat(mae,dilut,param_dict,dfchemo.index+273.15,
                                            sigma,growth_id,glc_up_id,prot_pool_id)
    solution["dataset"] = "chemostat"
    return  solution


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)
dfae_batch,dfan_batch = GEMS.load_exp_batch_data(os.path.join(path,'data/ExpGrowth.tsv'))
dfchemo = pd.read_csv(os.path.join(path,'data/Chemostat_exp_data.txt'),sep='\t',index_col=0)
sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]
Yobs_batch_an = {'data':dfan_batch.loc[sel_temp,'r_an'].values}



# In[ ]:


path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
simulated_results = pickle.load(open(file="../results/smcabc_gem_three_conditions_permuted_0_save_all_particles.pkl",mode='rb'))

# In[]


example_particle = simulated_results.all_particles[-1]




# In[]

example_model = pickle.load(open(os.path.join(path,'models/aerobic.pkl'),'rb'))



# %%
import timeit
start = timeit.timeit()
#fba_sol = GEMS.aerobic(example_particle)
fva_solution = aerobic_fva(example_particle)
end = timeit.timeit()
logging.info(end-start)
# %%
