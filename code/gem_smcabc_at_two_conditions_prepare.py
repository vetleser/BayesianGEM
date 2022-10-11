#!/usr/bin/env python
# coding: utf-8

# In[1]:


import abc_etc as abc
import numpy as np
import GEMS
import os
import pandas as pd
import pickle
import permute_parameters
import logging

# In[]

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s') 
# In[2]:


Yobs_batch = GEMS.aerobic_exp_data()
Yobs_chemo = GEMS.chemostat_exp_data()
#Yobs_batch_an = GEMS.anaerobic_exp_data()
dfae_batch,dfan_batch =GEMS.load_exp_batch_data('../data/ExpGrowth.tsv')
sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]
Yobs_batch_an = {'data':dfan_batch.loc[sel_temp,'r_an'].values}


# In[5]:


Yobs = {'rae':Yobs_batch['data'],
        'chemostat':Yobs_chemo['data'],
        'ran':Yobs_batch_an['data']}


# In[ ]:

overall_random_seed = 158
rng = np.random.default_rng(overall_random_seed)
path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)


candidate_frame: pd.DataFrame = pd.DataFrame(index=pd.DataFrame(range(4), names = ["simulation"])).reset_index()




# #### Define model settings

# In[ ]:


min_epsilon = -1.0 # equivalent to r2 score of 1
population_size = 100
outdir = '../results/permuted_smcabc_res'
if not os.path.exists(outdir):
    os.makedirs(outdir)
candidate_frame['outfile'] = [f'{outdir}/smcabc_gem_two_conditions_{simulation}_save_all_particles.pkl' for simulation in candidate_frame['simulation'])]
candidate_frame['random_seed'] = rng.choice(range(0,100000), candidate_frame.shape[0])

pickle.dump(file=open(file=f'{outdir}/simulation_skeleton.pkl',mode='wb'),obj=candidate_frame)

