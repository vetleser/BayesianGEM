# etcpy
#### Important: this is now only for ecYeast7 since the protein usage, NGAMT are different in other models. Will be adapted in the future if necessary.
This script provides several functions to map temperature effects on enzymes in
an enzyme constrained model. In the model, three temperature effects are considered:
* protein denaturation, which is described with a two-state denaturation model  
* enzyme kcat, which is described with a macromolecular theory macromolecular rate theory (MMRT): Hobbs J. et al., ACS Chem Biol. 2013  
* temperature dependent non-growth associated ATP maintenance (NGAM), which is from the experimentally measured NGAM at different temperatures under anaerobic conditions. The same dependence was assumed for aerobic growth  


## Usage:  
The input parameters:
* an enzyme constrained model in as `reframed.CBModel`
* a dictionary containing the following parameters of all enzymes in the model. Each reaction ID is a key in the dictionary. Each value in the dictionary is a dictionary having the following key-value pairs:
  * `dHTH`: the enthalpy change at convergence temperature TH (373.5 K), in J/mol
  * `dSTS`: the entropy change at convergence temperature TS (385 K), in J/mol/K
  * `dCpu`: the heat-capacity change of protein unfolding process, in J/mol/K
  * `dCpt`: the heat-capacity change between transition state and ground state in catalytic process, in J/mol/K
  * `Topt`: the optimal catalytic temperature of an enzymes, in K. At this temperature, the specific activity is maximal

```python
from etcpy import etc
```
1. temperature effect on protein denaturation  
```python
etc.reframed_mappers.map_fNT(model,T,param_dict,solver_instance=None)
```
model is the reframed model object. T is temperature in K. param_dict, the dictionary containing the following parameters of all enzymes in the model, solver_instance is an optional solver instance of type `reframed.solvers.Solver`. If a solver instance is provided, the solver instance will be updated in-place and the model argument will remain unchanged.

1. temperature effect on kcat  
```python
etc.reframed_mappers.map_kcatT(model,T,param_dict,solver_instance=None)
```
model is the reframed model object. T is temperature in K. param_dict, the dictionary containing the following parameters of all enzymes in the model, solver_instance is an optional solver instance of type `reframed.solvers.Solver`. If a solver instance is provided, the solver instance will be updated in-place and the model argument will remain unchanged.

3. temperature effect on NGAM
```python
etc.reframed_mappers.set_NGAMT(model,T)
```
model is the reframed model object **or** a solver instance


4. set sigma
```python
etc.reframed_mappers.set_sigma(model,sigma)
```

model is the reframed model object **or** a solver instance

5. simulate growth at different temperatures.
```python
etc.simulate_growth(model,Ts,sigma,param_dict)
```
Ts is a list of tempertures in K.

6. calculate dHTH, dSTS and dCpu of the denaturation process. There are two functions provided for two different scenarios.
```python 
dHTH,dSTS,dCpu = etc.thermal_parameters.get_dH_dS_dCpu_from_TmT90(Tm,T90)  #for a protein with experimental Tm and T90
dHTH,dSTS,dCpu = etc.thermal_parameters.get_dH_dS_dCpu_from_TmLength(Tm,proteinLength) # for protein with only Tm
```
Make sure that dCpu obtained from `etc.thermal_parameters.get_dH_dS_dCpu_from_TmT90(Tm,T90)` must be positive. If not, use `etc.thermal_parameters.get_dH_dS_dCpu_from_TmLength(Tm,proteinLength)`. 

1. Augment a dataframe `params` with the following columns: Tm,Tm_std,T90,dCpt,dCpt_std,Topt,Topt_std,Length with parameters from a dictionary `param_dict` containing some or all of the same parameters. This augmented dataframe is then converted to a dictionary of thermal parameters which can be used by `etc.simulate_growth`. Note that supplying an empty dictionary is equivalent to using the parameters in the dataframe as-is.
  ```python
  params = pd.read_csv('./model_enzyme_params.csv',index_col=0) # contains at least Tm,T90,dCpt,Topt,Length
  param_dict = etc.thermal_parameters.format_input(params, param_dict)
  ```


8. sample the uncertainties in the thermal parameters Tm,Topt and dCpt. Given a dataframe containing following columns:Tm,Tm_std,T90,dCpt,dCpt_std,Topt,Topt_std,Length. Randomly generate a new value from a normal distribution N(Tm,Tm_std) taking Tm as an example. 
```python 
new_params = etc.sample_data_uncertainty(params) # One can also specify columns to be sampled. The default is to sample all columns: [Tm,dCpt,Topt]
thermalparams = etc.calculate_thermal_params(new_params)
```
Then ```thermalparams``` could be used to simulate growth rate ```etc.simulate_growth(model,Ts,sigma,thermalparams)```


8. simulate chemostat data.  
(1) fix growth rate, set objective function as minimizing glucose uptatke rate  
(2) map temperature parameters  
(3) get minimal glucose uptake rate, then fix glucose uptake rate to this minimal value (\*1.001 for simulation purpose)  
(4) minimize enzyme usage
(5) return solution = `FBA(model)`
```python 
%%time
params = pd.read_csv('./model_enzyme_params.csv',index_col=0) # contains at least Tm,T90,dCpt,Topt,Length

param_dict = etc.thermal_parameters.format_input(params, {}) # Note the empty dictionary as the last argument
Ts = np.array([30,40,42,43,44,45,50])+273.15
growth_id = 'r_2111'
glc_up_id = 'r_1714_REV'
prot_pool_id = 'prot_pool_exchange'
solutions = etc.simulate_chemostat(model,0.1,param_dict,Ts,0.5,growth_id,glc_up_id,prot_pool_id)
```