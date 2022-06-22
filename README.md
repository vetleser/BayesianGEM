## Using Bayesian statistical learning to integrate temperature dependence in enzyme-constrained GEMs
<p align="center">
  <img  src="figures/logo.png">
</p>

#### Description of folders
* `code/` contains all scripts and detailed descrition can be found in `code/README.md`.
* `data/` contains all input data needed, including experimental and estimated thermal parameters.
* `models/` contains a list of yeast genome scale models with different settings used in this study.

#### Dependences
```
austin                      3.3.0
numpy                       1.21.0  
pandas                      1.3.3
scikit-learn                1.0.0
scipy                       1.7.1
reframed                    1.2.1
inspyred                    1.0.1
dill                        0.3.4  
jupyter                     1.0.0
matplotlib                  3.5.2
Gurobi                      9.1.2
```
The repository was tested with Python 3.8.12. The easiest way to install the dependencies is through the Conda package manager. Using Conda, the environment can be set up by running `conda env create --file .condaconfig.yml`.

#### Hardware
Since the Bayesian and evolutionary approach is computational expensive, all scripts except the visualization Jupyter notebook have to be done with a computer cluster. Those scripts have been designed for parallel computation through the use of SLURM. The visualization notebook take several seconds or minutes on a normal PC.

#### Reproduce the figures
(1) Clone this repository.  
(2) Install all required packages. This step takes at most several minutes.
(2) Download the pre-computed results from Figshare (https://zenodo.org/record/3996543#.X0J1BNP7S3I). Download the `results.tar.gz` file to the current directory and uncompress with 
```
tar -xzvf results.tar.gz
```
Then the figures in the manuscript can be reproduced by using Jupyter notebook

One can also recompute those results by following the introductions in `code/README.md` and the visualized by the above scripts.
