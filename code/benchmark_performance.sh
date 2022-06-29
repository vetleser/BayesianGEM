#!/bin/bash

# This script helps running benchmarks for particle evaluation for both reframed and cobrapy
# and uses git to which between the two frameworks

conda activate etcFBA
echo "Benchmarking with reframed"
python benchmark_performance.py
git checkout cfe99d7142e93ca12a36d8bb64d0a746d21c230f # Takes us back to the time the etcFBA framework used cobrapy
echo "Benchmarking with cobrapy"
python benchmark_performance.py
git switch -