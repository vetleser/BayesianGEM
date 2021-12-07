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

n_permutations = 3
random_seed = 158
rng = np.random.default_rng(random_seed)
path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)
permuted_params = [permute_parameters.permute_params(params, rng) for _ in range(n_permutations)]


# #### Define priors

# In[ ]:


priors = [dict() for _ in range(n_permutations)]

for i in range(n_permutations):
    for ind in params.index: 
        for col in ['Tm','Topt','dCpt']: 
            priors[i]['{0}_{1}'.format(ind,col)] = abc.RV('normal',
                                                          loc=permuted_params[i].loc[ind,col],
                                                          scale=permuted_params[i].loc[ind,col+'_std'])


# #### Define model settings

# In[ ]:


min_epsilon = -.9 # equivalent to r2 score of 1
population_size = 100
outfiles = [f'../results/smcabc_gem_three_conditions_permuted_{i}_save_all_particles.pkl' for i in range(n_permutations)]


# In[ ]:

models = []
for i, outfile in enumerate(outfiles):
    if not os.path.exists(outfile):
        logging.info('Initialize model')
        model = abc.SMCABC(GEMS.simulate_at_three_conditions_2,
                            priors[i],
                            min_epsilon,
                            population_size,
                            GEMS.distance_2,
                            Yobs,
                            outfile,
                            generation_size=128,
                            maxiter=100)
    else: model = pickle.load(open(outfile,'rb'))
    models.append(model)


# #### Run simulation

# In[ ]:

for i, model in enumerate(models):
    logging.info(f'Start simulations with prior parameter set {i}')
    model.run_simulation()

logging.info(f'DONE')
