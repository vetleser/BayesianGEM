import matplotlib
from matplotlib.font_manager import FontProperties
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
from matplotlib.lines import Line2D
import logging
import sa_etc_new as sa2

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')


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
    
def extract_simulated_annealing_final_generation_2(model: sa2.SimulatedAnnealing):
    return model.current_population

def contour_function(x,y):
    denomenator = 1+((x-1)**2+(y-1)**2)*((x+1)**2+(y+1)**2)
    R2 = 1/denomenator
    return -R2

def contour_function2(x,y):
    denomenator = 1+((x-1)**2+(y-1)**2)*((x+1)**2+(y+1)**2)*((x+1)**2+(y-1)**2)*((x-1)**2+(y+1)**2)
    R2 = 1/denomenator
    return -R2

def contour_function3(x, y):
    minima_list = [(-1, -1), (1, 1), (-1, 1), (1, -1)]
    sigma = 0.6
    R2 = np.zeros_like(x, dtype=float)

    for min_x, min_y in minima_list:
        dx = x - min_x
        dy = y - min_y
        exponent = -(dx**2 + dy**2) / (2 * sigma**2)
        R2 += np.exp(exponent)

    return -R2/1.0077467856174704

def contour_function4(x, y):
    minima_list = [(-1, 1), (1, 1), (-1, -1), (1, -1)]
    sigma = 0.6
    R2 = np.zeros_like(x, dtype=float)
    scaling_list = [1, 0.9, 0.8, 0.7]

    for i, (min_x, min_y) in enumerate(minima_list):
        dx = x - min_x
        dy = y - min_y
        exponent = -(dx**2 + dy**2) / (2 * sigma**2)
        R2 += np.exp(exponent) * scaling_list[i]

    return -R2/1.006582525974071
    
def contour_function5(x, y):
    minima_list = [(-1, 1), (1, 1), (-1, -1), (-3, 3)]
    sigma = 0.6
    R2 = np.zeros_like(x, dtype=float)
    scaling_list = [1, 0, 0, 1]

    for i, (min_x, min_y) in enumerate(minima_list):
        dx = x - min_x
        dy = y - min_y
        exponent = -(dx**2 + dy**2) / (2 * sigma**2)
        R2 += np.exp(exponent) * scaling_list[i]

    return -R2/1.0054182663306717

def general_contour_function(x, y, minima_list, scaling_list, sigma_list):
    logging.info(f"Minima list: {minima_list}, Scaling list: {scaling_list}, Sigma list: {sigma_list}")
    R2 = np.zeros_like(x, dtype=float)

    for i, (min_x, min_y) in enumerate(minima_list):
        dx = x - min_x
        dy = y - min_y
        exponent = -(dx**2 + dy**2) / (2 * sigma_list[i]**2)
        R2 += np.exp(exponent) * scaling_list[i]

    max_index = np.unravel_index(np.argmax(R2), R2.shape)
    max_x = x[max_index]
    max_y = y[max_index]
    max_val = R2[max_index]

    logging.info(f"Max R2 value: {max_val:.4f} at coordinates (x={max_x:.4f}, y={max_y:.4f})")

    max_R2 = np.max(R2)
    return -R2 / max_R2

def rastrigin_contour_function(x, y):
    A = 10
    n = 2
    # Calculate the Rastrigin function for x and y
    sum1 = (x**2 - A * np.cos(2 * np.pi * x))
    sum2 = (y**2 - A * np.cos(2 * np.pi * y))
    
    # Sum them together with the constant A * n
    return -(A * n + sum1 + sum2)/80.70658039



def make_plot(fig_type,
              final_temp_list, 
              minima_list,
              scaling_list,
              sigma_list,
              xlim=(-2.5, 2.5),
              ylim=(-2.5, 2.5),
              contour_func=None,
              filename = None,
            ):
    # Set up figure size and font
    fig = plt.figure(figsize=(15, 12.5))
    font = {'family': 'normal', 'weight': 'bold', 'size': 25}
    matplotlib.rc('font', **font)

    # Define the function to plot particles
    def plot_population_2(particles, simulation=0, cmap=None, markers=None):
        x, y = zip(*[(particle["x"], particle["y"]) for particle in particles])
        # Ensure unique color and marker per simulation
        plt.scatter(x, y, s=200, alpha=1, facecolors="none", edgecolors=cmap(simulation), 
                    marker=markers[simulation % len(markers)], linewidths=2.5)

        # CAN CHANGE THE LIMITS HERE
        # plt.xlim(xlim)
        # plt.ylim(ylim)
        plt.xlim(xlim)
        plt.ylim(ylim)

        plt.xlabel("x")
        plt.ylabel("y")

    # Number of simulations and plots (replace with actual values)
    n_plots = 1
    n_simulations = 4

    # Define the colormap and markers for simulations
    cmap = matplotlib.cm.get_cmap('tab10', n_simulations)
    markers = ["o", "s", "*", "^"]  # Customize with different markers if needed
    
    if contour_func is None:
        contour_func = lambda x, y: general_contour_function(x, y, minima_list, scaling_list, sigma_list)
    # Loop to create subplots for each plot
    for j, plot in enumerate(range(n_plots)):
        ax = plt.subplot(1, 1, j + 1)
        #X, Y = np.meshgrid(np.linspace(-10, 10, 1000), np.linspace(-10, 10, 1000))
        X, Y = np.meshgrid(np.linspace(xlim[0], xlim[1], 1000), np.linspace(ylim[0], ylim[1], 1000))
        
        Z = contour_func(X, Y) #Change contour function here
        ax.contour(X, Y, Z, 10, colors="black") 
        c = ax.pcolormesh(X, Y, Z, shading='auto', cmap='viridis', vmin=-1, vmax=0)
        cbar = fig.colorbar(c, ax=ax, label="Function value")
        cbar.set_ticks(np.arange(-1.0, 0.1, 0.2))  # optional: set custom ticks
        
        # Set subplot title
        ax.set_title(["A", "B", "C", "D"][j], loc="left", fontsize=45, fontweight="bold")

        # # Add text
        # plt.text(0.95, 0.95, f"Final Temperature: {[10, 1, 0.1, 0.01][j]}", fontsize=28, fontweight="normal", color="black", transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='right',     bbox=dict(
        #     facecolor='white',
        #     alpha=0.7,       # Transparency: 0 = fully transparent, 1 = opaque
        #     edgecolor='black',
        #     boxstyle='round,pad=0.3'
        # ))

            # Add text directly to the subplot
        ax.text(0.5, 1.0, f"Final Temperature: {final_temp_list[j]}",
                transform=ax.transAxes, fontsize=28, ha='center', va='bottom')


        # font_props = FontProperties(weight='normal', size=24)

        # # Add legend
        # plt.plot([], [], ' ', label=f"Final Temperature: {[10, 1, 0.1, 0.01][j]}")
        # plt.legend(loc='upper right', fontsize=24, frameon=True, prop=font_props)

        # Loop through all simulations for this plot and plot the population
        for i, simulation in enumerate(range(n_simulations)):
            # Get the final generation data for the given plot and simulation
            if filename == None:
                filename_sim = f"simanneal_toy_example_{plot}_{simulation}.pkl"
            else:
                filename_sim = f"{filename}_{plot}_{simulation}.pkl"
            model = load_pickle(f"../results/toy_example/{filename_sim}")
            particles = extract_simulated_annealing_final_generation(model)
            logging.info(f"Particles: {particles}")
            
            # Plot particles for this simulation
            plot_population_2(particles, simulation, cmap, markers)

    # Save and show the plot
    plt.tight_layout()
    logging.info(f"Saving figure")
    plt.savefig(f"../figures/toy_example/toy_example_{fig_type}2.png", dpi=300)
    plt.show()


logging.info("Gaussian fitness landscape, four minima, different scaling factors")
df = load_pickle("../results/toy_example/toy_example_df.pkl")
# make_plot(
#     fig_type="gaussian",
#     final_temp_list=[0.001],
#     minima_list=df["Minima_list"].tolist(),
#     scaling_list=df["Scaling_list"].tolist(),
#     sigma_list=df["Sigma_list"].tolist(),
#     xlim=(-5, 5),
#     ylim=(-5, 5),
#     filename="simanneal_toy_example",
# )

logging.info("Rastrigin fitness landscape")
make_plot(
    fig_type="rastr",
    final_temp_list=[0.001],
    minima_list=[],
    scaling_list=[],
    sigma_list=[],
    xlim=(-5.12, 5.12),
    ylim=(-5.12, 5.12),
    filename="simanneal_rastr",
    contour_func=rastrigin_contour_function
)

#logging.info("Four versions of simulated annealing, different final temperatures. Fitness function 3")
#Plotting!
# outdir = "../results/toy_example_new"
# n_simulations = 4
# n_plots = 4

# toy_example_df = load_pickle(f"{outdir}/toy_example_df.pkl")
# toy_example_df2 = load_pickle(f"{outdir}/toy_example_df2.pkl")

# make_plot(
#     fig_number="rastr",
#     final_temp_list=[1, 0.1, 0.01, 0.001],
#     minima_list=[(-1, 1), (1, 1), (-1, -1), (1, -1)],
#     scaling_list=[1, 0.9, 0.8, 0.7],
#     xlim=(-5.12, 5.12),
#     ylim=(-5.12, 5.12),
#     contour_func=rastrigin_contour_function,
#     filename="simanneal_rastr"
# )    

# for figure in toy_example_df["Figure"]:
#     row = toy_example_df.loc[toy_example_df["Figure"] == figure].iloc[0]  # robust row access
#     logging.info(f"Figure {figure} loaded")
#     make_plot(contour_func=rastrigin_contour_function,
#         fig_number=figure,
#         final_temp_list=row["Final_temp_list"],
#         minima_list=row["Minima_list"],
#         scaling_list=row["Scaling_list"],
#         xlim=row["X_lim"],
#         ylim=row["Y_lim"]
#     ) 
    #break #Uncomment this line to only run the first figure

# for figure in toy_example_df2["Figure"]:
#     row = toy_example_df2.loc[toy_example_df2["Figure"] == figure].iloc[0]  # robust row access
#     logging.info(f"Figure {figure} loaded")
#     make_plot(
#         fig_number=figure,
#         final_temp_list=row["Final_temp_list"],
#         minima_list=row["Minima_list"],
#         scaling_list=row["Scaling_list"],
#         xlim=row["X_lim"],
#         ylim=row["Y_lim"]
#     )

# toy_example_results = pd.DataFrame(index=pd.MultiIndex.from_product([range(n_plots),range(n_simulations)],names=["Plot","Simulation"])).reset_index()
# toy_example_results["modelfile"] = list(itertools.starmap(lambda plot,simulation: f"{'simanneal_test'}_{plot}_{simulation}.pkl",
#                                                        toy_example_results[["Plot","Simulation"]].itertuples(index=False,name=None))) #Can change name of file here
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

#     # CAN CHANGE THE LIMITS HERE
#     plt.xlim((-2.5, 2.5))
#     plt.ylim((-2.5, 2.5))
#     # plt.xlim((-5, 5))
#     # plt.ylim((-5, 5))
#     # plt.xlim((-5, 1))
#     # plt.ylim((-1, 5))

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
#     ax = plt.subplot(2, 2, j + 1)
#     X, Y = np.meshgrid(np.linspace(-10, 10, 1000), np.linspace(-10, 10, 1000))
#     Z = contour_function3(X, Y) #Change contour function here
#     ax.contour(X, Y, Z, 10, colors="black") 
    
#     # Set subplot title
#     ax.set_title(["A", "B", "C", "D"][j], loc="left", fontsize=45, fontweight="bold")

#     # # Add text
#     # plt.text(0.95, 0.95, f"Final Temperature: {[10, 1, 0.1, 0.01][j]}", fontsize=28, fontweight="normal", color="black", transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='right',     bbox=dict(
#     #     facecolor='white',
#     #     alpha=0.7,       # Transparency: 0 = fully transparent, 1 = opaque
#     #     edgecolor='black',
#     #     boxstyle='round,pad=0.3'
#     # ))

#         # Add text directly to the subplot
#     ax.text(0.5, 1.0, f"Final Temperature: {[1, 0.1, 0.01, 0.001][j]}",
#             transform=ax.transAxes, fontsize=28, ha='center', va='bottom')


#     # font_props = FontProperties(weight='normal', size=24)

#     # # Add legend
#     # plt.plot([], [], ' ', label=f"Final Temperature: {[10, 1, 0.1, 0.01][j]}")
#     # plt.legend(loc='upper right', fontsize=24, frameon=True, prop=font_props)

#     # Loop through all simulations for this plot and plot the population
#     for i, simulation in enumerate(range(n_simulations)):
#         # Get the final generation data for the given plot and simulation
#         particles = toy_example_results.loc[(plot, simulation), "final_generation"]
        
#         # Plot particles for this simulation
#         plot_population_2(particles, simulation, cmap, markers)

# # Save and show the plot
# plt.tight_layout()
# plt.savefig("../figures/toy_example/toy_example_test_wtitle.png", dpi=300)
# plt.show()

#----------------------------------------------------------------------------------

logging.info("Obtaining one Simulated Annealing and one Evolutionary Algorithm")

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


# #3 versions of simulated annealing, one version of evolutionary algorithm
# outdir = "../results/toy_example"
# n_simulations = 4
# n_plots = 4

# # toy_example_results = pd.DataFrame(index=pd.MultiIndex.from_product([range(n_plots),range(n_simulations)],names=["Plot","Simulation"])).reset_index()
# # toy_example_results["modelfile"] = list(itertools.starmap(lambda plot,simulation: f"{'simanneal'}_version_{plot}_plot_{simulation}.pkl",
# #                                                        toy_example_results[["Plot","Simulation"]].itertuples(index=False,name=None)))
# # toy_example_results["model"] = list(map(lambda filename: load_pickle(f"{outdir}/{filename}"),toy_example_results["modelfile"]))
# # toy_example_results["final_generation"] = list(itertools.starmap(lambda method, model: ( extract_simulated_annealing_final_generation)(model),
# #                                                toy_example_results[["Plot","model"]].itertuples(index=False,name=None)))
# # toy_example_results.set_index(["Plot","Simulation"],inplace=True)

# n_simulations = 4
# plots = [0, 1, 2, 3]

# # Create index for all (Plot, Simulation) combinations
# toy_example_results = pd.DataFrame(index=pd.MultiIndex.from_product([plots, range(n_simulations)], names=["Plot", "Simulation"])).reset_index()

# # Generate filenames
# def get_modelfile(plot, simulation):
#     if plot in [0, 1, 2]:
#         return f"simanneal_{plot}_{simulation}.pkl"
#     else:
#         return f"crowdingDE_{simulation}.pkl"

# toy_example_results["modelfile"] = list(itertools.starmap(get_modelfile, toy_example_results[["Plot", "Simulation"]].itertuples(index=False, name=None)))

# # Load models
# toy_example_results["model"] = list(map(lambda filename: load_pickle(f"{outdir}/{filename}"), toy_example_results["modelfile"]))

# # Extract final generation
# def get_final_generation(plot, model):
#     if plot in [0 ,1, 2]:
#         return extract_simulated_annealing_final_generation(model)
#     else:
#         return extract_evolution_final_generation(model)

# toy_example_results["final_generation"] = list(itertools.starmap(get_final_generation, toy_example_results[["Plot", "model"]].itertuples(index=False, name=None)))

# # Set index
# toy_example_results.set_index(["Plot", "Simulation"], inplace=True)
# print(toy_example_results.index)

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


# # Define the colormap and markers for simulations
# cmap = matplotlib.cm.get_cmap('tab10', n_simulations)
# markers = ["o", "s", "*", "^"]  # Customize with different markers if needed

# # Loop to create subplots for each plot
# for j, plot in enumerate(range(n_plots)):
#     plt.subplot(2, 2, j + 1)  # Create subplots in a 2x2 grid
#     X, Y = np.meshgrid(np.linspace(-5, 5, 1000), np.linspace(-5, 5, 1000))
#     Z = contour_function2(X, Y)
#     #levels = [-0.0001, -0.001, -0.01, -0.1].reverse()
#     plt.contour(X, Y, Z, levels=10, colors="black")  # Add contour lines
#     #plt.contourf(X, Y, Z, levels = 100, cmap = "viridis")  # Add contour lines
#     #plt.colorbar()
    

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
# plt.savefig("../figures/toy_example/toy_example_final_2.png", dpi=300)
# plt.show()


#----------------------------------------------------------------------------------
# #4 versions of simulated annealing, different min/max layers

# outdir = "../results/sa"
# n_simulations = 4
# n_plots = 4

# toy_example_results = pd.DataFrame(index=pd.MultiIndex.from_product([range(n_plots),range(n_simulations)],names=["Plot","Simulation"])).reset_index()
# toy_example_results["modelfile"] = list(itertools.starmap(lambda plot,simulation: f"{'simanneal_minmax_2'}_{plot}_{simulation}.pkl",
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
# # def plot_population_2(particles, simulation=0, cmap=None, markers=None):
# #     x, y = zip(*[(particle["x"], particle["y"]) for particle in particles])
# #     # Ensure unique color and marker per simulation
# #     plt.scatter(x, y, s=200, alpha=1, facecolors="none", edgecolors=cmap(simulation), 
# #                 marker=markers[simulation % len(markers)], linewidths=2.5)

# #     plt.xlim((-5, 5))
# #     plt.ylim((-5, 5))
# #     plt.xlabel("x")
# #     plt.ylabel("y")

# # Number of simulations and plots (replace with actual values)
# n_plots = 4
# n_simulations = 4

# # Define the colormap and markers for simulations
# cmap = matplotlib.cm.get_cmap('tab10', n_simulations)
# markers = ["o", "s", "*", "^"]  # Customize with different markers if needed

# # Loop to create subplots for each plot
# for j, plot in enumerate(range(n_plots)):
#     plt.subplot(2, 2, j + 1)  # Create subplots in a 2x2 grid
#     X, Y = np.meshgrid(np.linspace(-5, 5, 1000), np.linspace(-5, 5, 1000))
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
# plt.savefig("../figures/toy_example_simanneal_minmax_test.png", dpi=300)
# plt.show()
#----------------------------------------------------------------------------------

logging.info("Blank toy example fitness landscape")
df = load_pickle("../results/toy_example/toy_example_df.pkl")
# Set up figure size and font
plt.figure(figsize=(15, 12.5))
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
n_plots = 1
n_simulations = 4

# Define the colormap and markers for simulations
cmap = matplotlib.cm.get_cmap('tab10', n_simulations)
markers = ["o", "s", "*", "^"]  # Customize with different markers if needed

rastr_lim = (-5.12, 5.12)
plot_lim = (-10, 10)

# Loop to create subplots for each plot
#plt.subplot(2, 2, j + 1)  # Create subplots in a 2x2 grid
X, Y = np.meshgrid(np.linspace(plot_lim[0], plot_lim[1], 1000), np.linspace(plot_lim[0], plot_lim[1], 1000))
Z = rastrigin_contour_function(X, Y)
# Z = general_contour_function(X, Y,
#     minima_list=df["Minima_list"].tolist(),
#     scaling_list=df["Scaling_list"].tolist(),
#     sigma_list=df["Sigma_list"].tolist()
# )
plt.contour(X, Y, Z, 10, colors="black")  # Add contour lines
# plt.pcolormesh(X, Y, Z, shading='auto', cmap='viridis')
# plt.colorbar(label="Function value")
im = plt.imshow(Z, extent=(rastr_lim[0], rastr_lim[1], rastr_lim[0], rastr_lim[1]), origin='lower', cmap='viridis', aspect='auto')
# Add contour lines
cbar = plt.colorbar(im, label="Fitness Value")
cbar.set_ticks(np.arange(-1.0, 0.1, 0.2))  # ticks at 0.0, 0.2, 0.4, ..., 1.0




x_pos = [1, -1, 1, -1]
y_pos = [1, -1, -1, 1]


plt.title("Rastrigin Fitness Landscape", loc="center", fontsize=36, fontweight="bold")


plt.xlim(rastr_lim)
plt.ylim(rastr_lim)
plt.xlabel("x")
plt.ylabel("y")

# Save and show the plot
plt.tight_layout()
plt.savefig("../figures/toy_example/blank_rastr_3.png", dpi=300)
# plt.show()