#!/bin/bash

# This script helps running benchmarks for particle evaluation for both reframed and cobrapy
# and uses git to which between the two frameworks

conda activate etcFBA
echo "Benchmarking with reframed"
python benchmark_performance.py
git checkout cobrapy_profile # Takes us back to the time the etcFBA framework used cobrapy
echo "Benchmarking with COBRApy"
python benchmark_performance.py
git switch -
