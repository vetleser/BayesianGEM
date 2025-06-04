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

def general_contour_function(x, y, minima_list, scaling_list, sigma=0.6):
    R2 = np.zeros_like(x, dtype=float)

    for i, (min_x, min_y) in enumerate(minima_list):
        dx = x - min_x
        dy = y - min_y
        exponent = -(dx**2 + dy**2) / (2 * sigma**2)
        R2 += np.exp(exponent) * scaling_list[i]

    # Compute max R2 for normalization
    max_R2 = 0
    for i, (min_x, min_y) in enumerate(minima_list):
        exponent = -0 / (2 * sigma**2)  # dx=dy=0 at the minima point
        max_R2 += np.exp(exponent) * scaling_list[i]  # This sums up scaling_list[i]

    return -R2 / max_R2

def rastrigin_contour_function(x, y):
    A = 10
    n = 2
    # Calculate the Rastrigin function for x and y
    sum1 = (x**2 - A * np.cos(2 * np.pi * x))
    sum2 = (y**2 - A * np.cos(2 * np.pi * y))
    
    # Sum them together with the constant A * n
    return -(A * n + sum1 + sum2)/80.70658039



def make_plot(fig_number,
              final_temp_list, 
              minima_list,
              scaling_list,
              xlim=(-2.5, 2.5),
              ylim=(-2.5, 2.5),
              contour_func=None,
              filename = None,
            ):
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

        # CAN CHANGE THE LIMITS HERE
        # plt.xlim(xlim)
        # plt.ylim(ylim)
        plt.xlim(xlim)
        plt.ylim(ylim)

        plt.xlabel("x")
        plt.ylabel("y")

    # Number of simulations and plots (replace with actual values)
    n_plots = 4
    n_simulations = 4

    # Define the colormap and markers for simulations
    cmap = matplotlib.cm.get_cmap('tab10', n_simulations)
    markers = ["o", "s", "*", "^"]  # Customize with different markers if needed
    
    if contour_func is None:
        contour_func = lambda x, y: general_contour_function(x, y, minima_list, scaling_list)
    # Loop to create subplots for each plot
    for j, plot in enumerate(range(n_plots)):
        ax = plt.subplot(2, 2, j + 1)
        X, Y = np.meshgrid(np.linspace(-10, 10, 1000), np.linspace(-10, 10, 1000))
        
        Z = contour_func(X, Y) #Change contour function here
        ax.contour(X, Y, Z, 10, colors="black") 
        ax.pcolormesh(X, Y, Z, shading='auto', cmap='viridis')
        #plt.colorbar(label="Function value")
        
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
                filename_sim = f"simanneal_fig{fig_number}_new_{plot}_{simulation}.pkl"
            else:
                filename_sim = f"{filename}_{plot}_{simulation}.pkl"
            model = load_pickle(f"../results/toy_example_new/{filename_sim}")
            particles = extract_simulated_annealing_final_generation(model)
            logging.info(f"Particles: {particles}")
            
            # Plot particles for this simulation
            plot_population_2(particles, simulation, cmap, markers)

    # Save and show the plot
    plt.tight_layout()
    logging.info(f"Saving figure {fig_number}")
    plt.savefig(f"../figures/toy_example_new/toy_example_fig{fig_number}_new_wtitle.png", dpi=300)
    plt.show()




#logging.info("Four versions of simulated annealing, different final temperatures. Fitness function 3")
#Plotting!
outdir = "../results/toy_example_new"
n_simulations = 4
n_plots = 4

toy_example_df = load_pickle(f"{outdir}/toy_example_df.pkl")
toy_example_df2 = load_pickle(f"{outdir}/toy_example_df2.pkl")

make_plot(
    fig_number="rastr",
    final_temp_list=[1, 0.1, 0.01, 0.001],
    minima_list=[(-1, 1), (1, 1), (-1, -1), (1, -1)],
    scaling_list=[1, 0.9, 0.8, 0.7],
    xlim=(-5.12, 5.12),
    ylim=(-5.12, 5.12),
    contour_func=rastrigin_contour_function,
    filename="simanneal_rastr"
)    

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
# Set up figure size and font
plt.figure(figsize=(12.5, 12.5))
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

rastr_lim = (-5.12, 5.12)

# Loop to create subplots for each plot
#plt.subplot(2, 2, j + 1)  # Create subplots in a 2x2 grid
X, Y = np.meshgrid(np.linspace(rastr_lim[0], rastr_lim[1], 1000), np.linspace(rastr_lim[0], rastr_lim[1], 1000))
Z = rastrigin_contour_function(X, Y)
#plt.contour(X, Y, Z, 10, colors="black")  # Add contour lines
plt.pcolormesh(X, Y, Z, shading='auto', cmap='viridis')
plt.colorbar(label="Function value")



x_pos = [1, -1, 1, -1]
y_pos = [1, -1, -1, 1]
# Plot the four points
# for i in range(4):
#     plt.scatter(x_pos[i], y_pos[i], s=200, alpha=0.8, facecolors="green", edgecolors="green", 
#                 linewidths=2.5)
#     plt.text(x_pos[i], y_pos[i]+0.3, f"({x_pos[i]}, {y_pos[i]})", fontsize=22, fontweight="normal", color="black", 
#              verticalalignment='top', horizontalalignment='right', bbox=dict(
#         facecolor='white',
#         alpha=0.7,       # Transparency: 0 = fully transparent, 1 = opaque
#         edgecolor='black',
#         boxstyle='round,pad=0.3'
#     ))
# # plt.scatter(x_pos[i], y_pos[i], s=200, alpha=1, facecolors="none", edgecolors=cmap(i),
# plt.scatter(0, 0, s=200, alpha=0.8, facecolors="blue", edgecolors="blue")
# plt.text(0, 0+0.3, f"(-1.1)", fontsize=22, fontweight="normal", color="black", 
#          verticalalignment='top', horizontalalignment='right', bbox=dict(
#     facecolor='white',
#     alpha=0.7,       # Transparency: 0 = fully transparent, 1 = opaque
#     edgecolor='black',
#     boxstyle='round,pad=0.3'
# ))

# green_dot = Line2D([0], [0], marker='o', color='w', label='Global minima, R² = -1.0',
#                        markerfacecolor='green', markersize=15)
# blue_dot = Line2D([0], [0], marker='o', color='w', label='Local minimum, R² ≈ -0.25',
#                       markerfacecolor='blue', markersize=15)

# # Add legend
# plt.legend(handles=[green_dot, blue_dot], loc='upper right', fontsize=14, frameon=True)

plt.title("Toy Example Rastrigin Fitness Landscape", loc="center", fontsize=36, fontweight="bold")


plt.xlim(rastr_lim)
plt.ylim(rastr_lim)
plt.xlabel("x")
plt.ylabel("y")

# Save and show the plot
plt.tight_layout()
plt.savefig("../figures/toy_example_new/toy_example_blank_rastrigin.png", dpi=300)
plt.show()