import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.cluster as cluster
import pickle
import GEMS
import itertools
import dill
import abc_etc
import evo_etc
import sa_etc

# Convenient pickle wrappers
def load_pickle(filename):
    return dill.load(open(file=filename,mode='rb'))
        
def dump_pickle(obj,filename):
    return dill.dump(obj=obj,file=open(file=filename, mode='wb'))

def extract_bayesian_final_generation(model: abc_etc.SMCABC ,population_size=128):
    particles_to_pick = np.argsort(model.all_distances)[:population_size]
    return [model.all_particles[particle] for particle in particles_to_pick]

def extract_evolution_final_generation(model: evo_etc.CrowdingDE):
    particles_to_pick = model.population[-1]
    return [model.all_particles[particle] for particle in particles_to_pick]

def extract_simulated_annealing_final_generation(model: sa_etc.SimulatedAnnealing):
    try:
        particles_to_pick = model.population[-1]
        return [model.all_particles[particle] for particle in particles_to_pick]
    except IndexError:
        return model.population_old

def contour_function(x,y):
    denomenator = 1+((x-1)**2+(y-1)**2)*((x+1)**2+(y+1)**2)
    R2 = 1/denomenator
    return -R2
    

#Obtaining one Simulated Annealing and one Evolutionary Algorithm
# outdir = "../results/sa"
# n_simulations = 4
# toy_example_results = pd.DataFrame(index=pd.MultiIndex.from_product([["Simulated Annealing","Evolutionary"],range(n_simulations)],names=["Method","Simulation"])).reset_index()
# toy_example_results["modelfile"] = list(itertools.starmap(lambda method,simulation: f"{'simanneal' if method == 'Simulated Annealing' else 'crowdingDE'}_{simulation}.pkl",
#                                                       toy_example_results[["Method","Simulation"]].itertuples(index=False,name=None)))
# toy_example_results["model"] = list(map(lambda filename: load_pickle(f"{outdir}/{filename}"),toy_example_results["modelfile"]))
# toy_example_results["final_generation"] = list(itertools.starmap(lambda method, model: (extract_simulated_annealing_final_generation if method == "Simulated Annealing" else extract_evolution_final_generation)(model),
#                                               toy_example_results[["Method","model"]].itertuples(index=False,name=None)))
# toy_example_results.set_index(["Method","Simulation"],inplace=True)




# plt.figure(figsize=(2*12.5,12.5))
# font = {'family' : 'normal',
#         'weight' : 'bold',
#         'size'   : 25}
# matplotlib.rc('font', **font)
# def plot_population(particles,simulation=0):
#     cmap = matplotlib.cm.get_cmap('tab10',n_simulations)
#     markers = ["o","s","*", "^"]
#     x,y = zip(*[(particle["x"],particle["y"]) for particle in particles])
#     plt.scatter(x,y,s=200,alpha=1,facecolors="none",edgecolors=cmap(simulation),marker=markers[simulation % 4],linewidths=2.5)
#     plt.xlim((-2,2))
#     plt.ylim((-2,2))
#     plt.xlabel("x")
#     plt.ylabel("y")


# for i, method in enumerate(["Simulated Annealing","Evolutionary"],start=0):
#     plt.subplot(1,2,i+1)
#     X, Y = np.meshgrid(np.linspace(-2,2,1000),np.linspace(-2,2,1000))
#     Z = contour_function(X,Y)
#     plt.contour(X,Y,Z,10,colors="black")
#     plt.title(["A","B"][i],loc="left",fontsize=45,fontweight="bold")
#     for j, simulation in enumerate(toy_example_results.loc[method].index):
#         plot_population(toy_example_results.loc[(method,simulation),"final_generation"],j)
# plt.savefig("../figures/toy_example_3.png",dpi=300)
# plt.show()

#--------------------------------------------------------------------------------------------




#Plotting!
# outdir = "../results/sa"
# n_simulations = 4
# n_plots = 4

# toy_example_results = pd.DataFrame(index=pd.MultiIndex.from_product([range(n_plots),range(n_simulations)],names=["Plot","Simulation"])).reset_index()
# toy_example_results["modelfile"] = list(itertools.starmap(lambda plot,simulation: f"{'simanneal'}_{plot}_{simulation}.pkl",
#                                                        toy_example_results[["Plot","Simulation"]].itertuples(index=False,name=None)))
# toy_example_results["model"] = list(map(lambda filename: load_pickle(f"{outdir}/{filename}"),toy_example_results["modelfile"]))
# toy_example_results["final_generation"] = list(itertools.starmap(lambda method, model: ( extract_simulated_annealing_final_generation)(model),
#                                                toy_example_results[["Plot","model"]].itertuples(index=False,name=None)))
# toy_example_results.set_index(["Plot","Simulation"],inplace=True)



# # Set up figure size and font
# plt.figure(figsize=(2 * 12.5, 2 * 12.5))
# font = {'family': 'normal', 'weight': 'bold', 'size': 25}
# matplotlib.rc('font', **font)

# # Define the function to plot particles
# def plot_population_2(particles, simulation=0, cmap=None, markers=None):
#     x, y = zip(*[(particle["x"], particle["y"]) for particle in particles])
#     # Ensure unique color and marker per simulation
#     plt.scatter(x, y, s=200, alpha=1, facecolors="none", edgecolors=cmap(simulation), 
#                 marker=markers[simulation % len(markers)], linewidths=2.5)

#     plt.xlim((-2, 2))
#     plt.ylim((-2, 2))
#     plt.xlabel("x")
#     plt.ylabel("y")

# # Number of simulations and plots (replace with actual values)
# n_plots = 4
# n_simulations = 4

# # Define the colormap and markers for simulations
# cmap = matplotlib.cm.get_cmap('tab10', n_simulations)
# markers = ["o", "s", "*", "^"]  # Customize with different markers if needed

# # Loop to create subplots for each plot
# for j, plot in enumerate(range(n_plots)):
#     plt.subplot(2, 2, j + 1)  # Create subplots in a 2x2 grid
#     X, Y = np.meshgrid(np.linspace(-2, 2, 1000), np.linspace(-2, 2, 1000))
#     Z = contour_function(X, Y)
#     plt.contour(X, Y, Z, 10, colors="black")  # Add contour lines
    
#     # Set subplot title
#     plt.title(["A", "B", "C", "D"][j], loc="left", fontsize=45, fontweight="bold")
    
#     # Loop through all simulations for this plot and plot the population
#     for i, simulation in enumerate(range(n_simulations)):
#         # Get the final generation data for the given plot and simulation
#         particles = toy_example_results.loc[(plot, simulation), "final_generation"]
        
#         # Plot particles for this simulation
#         plot_population_2(particles, simulation, cmap, markers)

# # Save and show the plot
# plt.tight_layout()
# plt.savefig("../figures/toy_example_four_simanneal_new.png", dpi=300)
# plt.show()

#----------------------------------------------------------------------------------

#3 versions of simulated annealing, one version of evolutionary algorithm
outdir = "../results/sa"
n_simulations = 4
n_plots = 4

# toy_example_results = pd.DataFrame(index=pd.MultiIndex.from_product([range(n_plots),range(n_simulations)],names=["Plot","Simulation"])).reset_index()
# toy_example_results["modelfile"] = list(itertools.starmap(lambda plot,simulation: f"{'simanneal'}_version_{plot}_plot_{simulation}.pkl",
#                                                        toy_example_results[["Plot","Simulation"]].itertuples(index=False,name=None)))
# toy_example_results["model"] = list(map(lambda filename: load_pickle(f"{outdir}/{filename}"),toy_example_results["modelfile"]))
# toy_example_results["final_generation"] = list(itertools.starmap(lambda method, model: ( extract_simulated_annealing_final_generation)(model),
#                                                toy_example_results[["Plot","model"]].itertuples(index=False,name=None)))
# toy_example_results.set_index(["Plot","Simulation"],inplace=True)

n_simulations = 4
plots = [1, 2, 3, 4]

# Create index for all (Plot, Simulation) combinations
toy_example_results = pd.DataFrame(index=pd.MultiIndex.from_product([plots, range(n_simulations)], names=["Plot", "Simulation"])).reset_index()

# Generate filenames
def get_modelfile(plot, simulation):
    if plot in [1, 2, 3]:
        return f"simanneal_version_{plot}_plot_{simulation}.pkl"
    else:
        return f"crowdingDE_{simulation}.pkl"

toy_example_results["modelfile"] = list(itertools.starmap(get_modelfile, toy_example_results[["Plot", "Simulation"]].itertuples(index=False, name=None)))

# Load models
toy_example_results["model"] = list(map(lambda filename: load_pickle(f"{outdir}/{filename}"), toy_example_results["modelfile"]))

# Extract final generation
def get_final_generation(plot, model):
    if plot in [1, 2, 3]:
        return extract_simulated_annealing_final_generation(model)
    else:
        return extract_evolution_final_generation(model)

toy_example_results["final_generation"] = list(itertools.starmap(get_final_generation, toy_example_results[["Plot", "model"]].itertuples(index=False, name=None)))

# Set index
toy_example_results.set_index(["Plot", "Simulation"], inplace=True)
print(toy_example_results.index)

# Set up figure size and font
plt.figure(figsize=(2 * 12.5, 2 * 12.5))
font = {'family': 'normal', 'weight': 'bold', 'size': 25}
matplotlib.rc('font', **font)

# Define the function to plot particles
def plot_population_2(particles, simulation=0, cmap=None, markers=None):
    x, y = zip(*[(particle["x"], particle["y"]) for particle in particles])
    # Ensure unique color and marker per simulation
    plt.scatter(x, y, s=200, alpha=1, facecolors="none", edgecolors=cmap(simulation), 
                marker=markers[simulation % len(markers)], linewidths=2.5)

    plt.xlim((-5, 5))
    plt.ylim((-5, 5))
    plt.xlabel("x")
    plt.ylabel("y")


# Define the colormap and markers for simulations
cmap = matplotlib.cm.get_cmap('tab10', n_simulations)
markers = ["o", "s", "*", "^"]  # Customize with different markers if needed

# Loop to create subplots for each plot
for j, plot in enumerate(range(n_plots)):
    plt.subplot(2, 2, j + 1)  # Create subplots in a 2x2 grid
    X, Y = np.meshgrid(np.linspace(-5, 5, 1000), np.linspace(-5, 5, 1000))
    Z = contour_function(X, Y)
    plt.contour(X, Y, Z, 10, colors="black")  # Add contour lines
    
    # Set subplot title
    plt.title(["A", "B", "C", "D"][j], loc="left", fontsize=45, fontweight="bold")
    
    # Loop through all simulations for this plot and plot the population
    for i, simulation in enumerate(range(n_simulations)):
        # Get the final generation data for the given plot and simulation
        particles = toy_example_results.loc[(plot+1, simulation), "final_generation"]
        
        # Plot particles for this simulation
        plot_population_2(particles, simulation, cmap, markers)

# Save and show the plot
plt.tight_layout()
plt.savefig("../figures/toy_example_final_2.png", dpi=300)
plt.show()


#----------------------------------------------------------------------------------
#4 versions of simulated annealing, different min/max layers

outdir = "../results/sa"
n_simulations = 4
n_plots = 4

toy_example_results = pd.DataFrame(index=pd.MultiIndex.from_product([range(n_plots),range(n_simulations)],names=["Plot","Simulation"])).reset_index()
toy_example_results["modelfile"] = list(itertools.starmap(lambda plot,simulation: f"{'simanneal_minmax_2'}_{plot}_{simulation}.pkl",
                                                       toy_example_results[["Plot","Simulation"]].itertuples(index=False,name=None)))
toy_example_results["model"] = list(map(lambda filename: load_pickle(f"{outdir}/{filename}"),toy_example_results["modelfile"]))
toy_example_results["final_generation"] = list(itertools.starmap(lambda method, model: ( extract_simulated_annealing_final_generation)(model),
                                               toy_example_results[["Plot","model"]].itertuples(index=False,name=None)))
toy_example_results.set_index(["Plot","Simulation"],inplace=True)

# Set up figure size and font
plt.figure(figsize=(2 * 12.5, 2 * 12.5))
font = {'family': 'normal', 'weight': 'bold', 'size': 25}
matplotlib.rc('font', **font)

# Define the function to plot particles
# def plot_population_2(particles, simulation=0, cmap=None, markers=None):
#     x, y = zip(*[(particle["x"], particle["y"]) for particle in particles])
#     # Ensure unique color and marker per simulation
#     plt.scatter(x, y, s=200, alpha=1, facecolors="none", edgecolors=cmap(simulation), 
#                 marker=markers[simulation % len(markers)], linewidths=2.5)

#     plt.xlim((-5, 5))
#     plt.ylim((-5, 5))
#     plt.xlabel("x")
#     plt.ylabel("y")

# Number of simulations and plots (replace with actual values)
n_plots = 4
n_simulations = 4

# Define the colormap and markers for simulations
cmap = matplotlib.cm.get_cmap('tab10', n_simulations)
markers = ["o", "s", "*", "^"]  # Customize with different markers if needed

# Loop to create subplots for each plot
for j, plot in enumerate(range(n_plots)):
    plt.subplot(2, 2, j + 1)  # Create subplots in a 2x2 grid
    X, Y = np.meshgrid(np.linspace(-5, 5, 1000), np.linspace(-5, 5, 1000))
    Z = contour_function(X, Y)
    plt.contour(X, Y, Z, 10, colors="black")  # Add contour lines
    
    # Set subplot title
    plt.title(["A", "B", "C", "D"][j], loc="left", fontsize=45, fontweight="bold")
    
    # Loop through all simulations for this plot and plot the population
    for i, simulation in enumerate(range(n_simulations)):
        # Get the final generation data for the given plot and simulation
        particles = toy_example_results.loc[(plot, simulation), "final_generation"]
        
        # Plot particles for this simulation
        plot_population_2(particles, simulation, cmap, markers)

# Save and show the plot
plt.tight_layout()
plt.savefig("../figures/toy_example_simanneal_minmax_test.png", dpi=300)
plt.show()
#----------------------------------------------------------------------------------