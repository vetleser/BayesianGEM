#!/usr/bin/env python
# coding: utf-8

# In[1]:


from random_sampler import RV
import numpy as np
import GEMS
import os
import pandas as pd
import pickle
import logging

# In[]

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s') 

# In[ ]:


final_temp = [0.0001]
#normalize = [False]
temp_step_size = [0.1, 0.5]
dCpt_step_size = [0.5, 1.0]
move_type = ['normal']

n_replicates = 2
overall_random_seed = 200
rng = np.random.default_rng(overall_random_seed)
path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)

candidate_frame = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [final_temp, move_type, temp_step_size, dCpt_step_size, range(n_replicates)],
        names=["final_temp", "move_type", "temp_step_size", "dCpt_step_size", "simulation"]
    )
).reset_index()



# In[ ]:


min_epsilon = -1.0 # equivalent to r2 score of 1
population_size = 100
outdir = '../results/sa'
if not os.path.exists(outdir):
    os.makedirs(outdir)
candidate_frame['outfile'] = [f'{outdir}/smcsa_gem_june2_{temp_step_size}_{dCpt_step_size}_{simulation}.pkl' for temp_step_size,
 dCpt_step_size, simulation in zip(candidate_frame['temp_step_size'], candidate_frame['dCpt_step_size'],candidate_frame['simulation'])]
candidate_frame['random_seed'] = rng.choice(range(0,100000), candidate_frame.shape[0])

pickle.dump(file=open(file=f'{outdir}/simulation_skeleton.pkl',mode='wb'),obj=candidate_frame)

logging.info(f'Frame columns : {candidate_frame.columns}')
logging.info(f'Simulation skeleton : {candidate_frame.to_string()}')

