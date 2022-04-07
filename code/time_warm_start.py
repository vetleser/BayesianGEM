#!/usr/bin/env python
# coding: utf-8
# The purpose of this script is to benchmark warm-start FBA versus
# the naive approach of optimizing any individual subproblem
# 
# In[]
import time
from timeit import timeit
import GEMS
from etcpy import etc
import os
import pandas as pd
import logging
import abc_etc as abc
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

path = '../'
params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)
particle = dict()
for ind in params.index: 
    for col in ['Tm','Topt','dCpt']: 
        particle['{0}_{1}'.format(ind,col)] = params.loc[ind,col]

# In[]
# tic = time.perf_counter()
# aerobic_cold_start = GEMS.aerobic(particle,warm_start=False)
# toc = time.perf_counter()
# print(f"Computed aerobic conditions without warm start in {toc-tic} seconds.")
tic = time.perf_counter()
total_cold_start = GEMS.simulate_at_three_conditions_2(particle,warm_start=False)
toc = time.perf_counter()
print(f"Computed results at three conditions without warm start in {toc-tic} seconds.")

# In[]
def junk():
    # In[]
    tic = time.perf_counter()
    aerobic_cold_start = GEMS.aerobic(particle,warm_start=False)
    toc = time.perf_counter()
    print(f"Computed aerobic conditions without warm start in {toc-tic} seconds.")


    # In[]
    tic = time.perf_counter()
    aerobic_warm_start = GEMS.aerobic(particle,warm_start=True)
    toc = time.perf_counter()
    print(f"Computed aerobic conditions with warm start in {toc-tic} seconds.")

    # In[]
    tic = time.perf_counter()
    total_cold_start = GEMS.simulate_at_three_conditions_2(particle,warm_start=False)
    toc = time.perf_counter()
    print(f"Computed results at three conditions without warm start in {toc-tic} seconds.")

    # In[]
    tic = time.perf_counter()
    total_warm_start = GEMS.simulate_at_three_conditions_2(particle,warm_start=True)
    toc = time.perf_counter()
    print(f"Computed results at three conditions with warm start in {toc-tic} seconds.")
    # In[]
    tic = time.perf_counter()
    chemostat_warm_start = GEMS.chemostat(particle,warm_start=True)
    toc = time.perf_counter()
    print(f"Computed chemostat results with warm start in {toc-tic} seconds.")