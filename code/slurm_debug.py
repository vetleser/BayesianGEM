#!/usr/bin/env python
# coding: utf-8

import pickle
import os

task_idx = int(os.environ["SLURM_ARRAY_TASK_ID"])

pickle.dump(task_idx,file=open(f"debug_results_{task_idx}.pkl",'wb'))
