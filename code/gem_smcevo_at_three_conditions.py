#!/usr/bin/env python
# coding: utf-8

# In[1]:


import dill
import evo_etc as evo
import numpy as np
import GEMS
import os
import pandas as pd
import pickle
import logging

# In[]

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s') 

# In[2]:


def main():
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


    path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
    params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)


    # #### Define priors

    # In[ ]:


    priors = dict()
    for ind in params.index: 
        for col in ['Tm','Topt','dCpt']: 
            priors['{0}_{1}'.format(ind,col)] = evo.RV('normal',
                                                        loc=params.loc[ind,col],
                                                        scale=params.loc[ind,col+'_std'])


    # #### Define model settings

    # In[ ]:

    rng = evo.default_rng(39322)
    min_epsilon = -.9 # equivalent to r2 score of 1
    population_size = 100
    outfile = '../results/smcevo_gem_three_conditions_save_all_particles.pkl'


    # In[ ]:


    if not os.path.exists(outfile):
        logging.info('Initialize model')
        model = evo.GA(simulator= GEMS.simulate_at_three_conditions_2,
                            priors=priors,
                            min_epsilon=min_epsilon,
                            generation_size=population_size,
                            distance_function=GEMS.distance_2,
                            Yobs=Yobs,
                            outfile=outfile,
                            maxiter=400,
                            rng=rng,
                            mutation_frequency=20,
                            mutation_prob=.5)
    else: 
        model = dill.load(open(outfile,'rb'))


    # #### Run simulation

    # In[ ]:


    logging.info('start simulations')
    model.run_simulation()


if __name__ == '__main__':
    main()
