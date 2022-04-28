#!/usr/bin/env python
# coding: utf-8

# The script is supposed to finish a job which randomly crashed in the 491th iteration
# The exact random state of the last iteration cannot be exactly reproduced, 
# but the real impact should be minimal

#In[]
import abc_etc as abc
import numpy as np
import GEMS
import os
import pandas as pd
import pickle
import logging

# In[]

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
task_idx = 1 # Job at array index is to be resumed




# In[5]:



outdir = '../results/permuted_smcabc_res'
candidate_frame: pd.DataFrame = pickle.load(file=open(file=f'{outdir}/simulation_skeleton.pkl',mode='rb'))
entry = candidate_frame.iloc[task_idx]



# In[]

min_epsilon = -1.0 # equivalent to r2 score of 1
population_size = 100
logging.info('Initialize model')


priors, origin, status, outfile, random_seed = entry[["priors", "origin", "status", "outfile","random_seed"]]
model = pickle.load(open(f"{outdir}/smcabc_gem_three_conditions_updated_permuted_0_replicate_save_all_particles.pkl",'rb'))



# In[ ]:
# #### Run simulation

# np.random.seed(random_seed)
logging.info(f'Resuming simulations with prior parameter set {origin}, {status}')
model.run_simulation()

logging.info(f'DONE')
