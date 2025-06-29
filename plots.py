"""
In this file, various plotting helper functions are defined.
"""

import matplotlib.pyplot as plt
import numpy as np
import units_conversion

########################################################   
############### Energies  ##############################
########################################################

def plot_energies(number_time_steps, time_step, kinetic_energy, potential_energy, total_energy=None):
    """
    This functions creates a plot of the evolution of the total energy, kinetic energy and potential energy of the system
    
    Parameters
    ----------
    number_time_steps: int
        The number of time steps of the simulation
    time_step: float
        The time step used in the simulation
    kinetic_energy: np.array
        Contains kinetic energies over time
        Array dimensions: number_time_steps
    potential_energy: np.array
        Contains potential energies over time
        Array dimensions: number_time_steps
    total_energy: np.array 
        Contains total energies over time
        Array dimensions: number_time_stepsof number_time_steps elements
        
    Returns 
    ----------
    """

    # Compute total energy if not provided
    if total_energy is None:
        total_energy = kinetic_energy + potential_energy

    # Define an array of times of the same dimension as the energy arrays
    times = np.arange(0, number_time_steps*time_step, time_step)

    # Plot energies
    plt.figure(figsize=(8, 5))
    plt.plot(times, kinetic_energy, label="Kinetic Energy", linestyle="--")
    plt.plot(times, potential_energy, label="Potential Energy", linestyle=":")
    plt.plot(times, total_energy, label="Total Energy", linestyle="-")
    
    plt.xlabel("Time")
    plt.ylabel("Energy")
    plt.legend()
    plt.title("Energy Evolution Over Time")
    plt.grid()
    plt.show()
    
########################################################   
############### Velocities #############################
########################################################

def velocity_distribution(velocities, temp):
    """
    This functions creates a plot of the inital velocity distribution in each coordinate
    
    Parameters
    ----------
    velocities: np.array
        Contains initial velocities
        array dimension: number_atoms x 3
    temp: float
        Desired temperature according to which Boltzmann distribution is initialized (for reference plotting)
        
    Returns 
    ----------
    """

    # Plot initial velocity distributions
    fig, ax = plt.subplots(1,3)
    
    # create a probability density in order to assure commensurability with normal distribution
    distance_histogram_x, bins_x = np.histogram(velocities[:,0], bins=30, density=True)        
    distance_histogram_y, bins_y = np.histogram(velocities[:,1], bins=30, density=True)
    distance_histogram_z, bins_z = np.histogram(velocities[:,2], bins=30, density=True)

    # determine value of Gaussian distribution at midpoints
    gaussian_x = 1/np.sqrt(2*np.pi*temp)*np.exp(-(((bins_x[:-1]+bins_x[1:])/2)**2)/(2*temp))  
    gaussian_y = 1/np.sqrt(2*np.pi*temp)*np.exp(-(((bins_y[:-1]+bins_x[1:])/2)**2)/(2*temp))  
    gaussian_z = 1/np.sqrt(2*np.pi*temp)*np.exp(-(((bins_z[:-1]+bins_x[1:])/2)**2)/(2*temp))  

    # plot initializations and comparison with Gaussian distribution
    ax[0].plot(bins_x[:-1], distance_histogram_x)
    ax[0].plot(bins_x[:-1], gaussian_x)
    ax[1].plot(bins_y[:-1], distance_histogram_y)
    ax[1].plot(bins_y[:-1], gaussian_y)
    ax[2].plot(bins_z[:-1], distance_histogram_z)
    ax[2].plot(bins_z[:-1], gaussian_z)

    ax[0].set_xlabel("x-velocities")
    ax[1].set_xlabel("y-velocities")
    ax[2].set_xlabel("z-velocities")

    ax[0].set_ylabel("Probability")
    ax[0].set_title("x-vel initialization")
    ax[1].set_title("y-vel initialization")
    ax[2].set_title("z-vel initialization")

    fig.suptitle("Velocity Component Distributions vs. Gaussian", fontsize=16)

    plt.show()
    
########################################################   
####### Intergration method comparison #################
########################################################

def euler_vs_verlet(number_time_steps, time_step, total_energy_euler, total_energy_verlet):
    """
    This functions creates a plot of the evolution of the total energy computed with the Euler and Velocity-Verlet methods
    
    Parameters
    ----------    
    number_time_steps: int
        The number of time steps of the simulation
    time_step: float
        The time step used in the simulation
    total_energy_euler: np.array
        Contains total energies according to Euler integration
        Array dimension: number_time_steps
    total_energy_verlet: np.array
        Contains total energies according to Velocity-Verlet integration
        Array dimension: number_time_steps
        
    Returns  
    ----------   
    """

    # Define an array of times of the same dimension as the energy arrays
    times = np.arange(0, number_time_steps*time_step, time_step)

    # Plot energies
    plt.figure(figsize=(8, 5))
    plt.plot(times, total_energy_euler, label="Total Energy (Euler)", linestyle="-")
    plt.plot(times, total_energy_verlet, label="Total Energy (Verlet)", linestyle="--")
    
    plt.xlabel("Time")
    plt.ylabel("Energy")
    plt.legend()
    plt.title("Comparison of energy conservation using Euler and Velocity-Verlet methods")
    plt.grid()
    plt.show()

########################################################   
################# Temperature  #########################
########################################################

def plot_temp(number_time_steps, time_step, temperature, normal_temp, SI):
    """
    This functions creates a plot of the evolution of the temperature of the system
    
    Parameters
    ----------    
    number_time_steps: int
        The number of time steps of the simulation
    time_step: float
        The time step used in the simulation
    temperature: np.array
        Contains temperatures at each timestep
        Array dimension: number_time_steps
    normal_temp: float
        Desired system temperature
    SI: Boolean
        Specifies if the plot is using normalised or SI temperatures

    Returns  
    ----------   
    """
    
    if SI==False:
        # Define an array of times of the same dimension as the energy arrays
        times = np.arange(0, number_time_steps*time_step, time_step)
    
        # Plot energies
        plt.figure(figsize=(8, 5))
        plt.plot(times, temperature)
        plt.xlabel("Time")

        plt.ylabel("Normalised Temperature")

    else:
        # Define an array of times of the same dimension as the energy arrays
        times = np.arange(0, number_time_steps*time_step, time_step)
        times=units_conversion.dimensionless_to_SI(times, "time")
        
        # Plot energies
        plt.figure(figsize=(8, 5))
        plt.plot(times, temperature)
        plt.xlabel("Time (s)")
        
        plt.ylabel("Temperature (K)")

    plt.axhline(normal_temp, color = "green", linestyle="--")
    plt.title("Temperature Evolution Over Time")
    plt.grid()
    plt.show()
    
    
########################################################   
################# Pair Correlation #####################
########################################################

def plot_pair_correlation(bins, function_vals):
    """
    This functions creates a plot of the pair correlation, from which the phase of 
    the system can be deduced.
    
    Parameters
    ----------    
    bins: int
        The number of bins for distances
    function_vals: np.array
        Pair correlation values that indicate difference to ideal gas
        Array dimension: number of bins
        
    Returns  
    ----------   
    """
    
    # Plot pair correlation function
    plt.figure(figsize=(8, 5))
    plt.plot(bins, function_vals)
    plt.xlabel("r", fontsize = 16)
    plt.ylabel("g(r)", fontsize = 16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.axhline(1, color = "green", linestyle="--")
    plt.title("Radial distribution function g(r) in Argon MD simulation", fontsize = 18)
    plt.grid()
    plt.show()

########################################################   
################# Diffusion ############################
########################################################

def plot_diffusion(time_step, mean_squared_displacement, diffusion):
    """
    This functions creates a plot of the diffusion.
    
    Parameters
    ----------    

    time_step: float
        The time step used in the simulation
    mean_squared_displacement: np.array
        Displacement values over time
        Array dimension: number_of_time_steps
    diffusion: np.array
        Diffusion values over time
        Array dimension: number_of_time_steps
        
    Returns  
    ----------   
    """
    
    # Construct time vector starting from initial equilibrium time
    number_time_steps = mean_squared_displacement.shape[0]
    times = np.arange(0, number_time_steps*time_step, time_step)
    
    # Plot the mean_squared_displacement
    plt.figure(figsize=(8,5))
    plt.plot(times, mean_squared_displacement)
    plt.xlabel("Time")
    plt.ylabel("Mean squared displacement")
    plt.title("Evolution of the mean squared displacement")
    plt.grid()
    plt.show()

    # Plot the diffusion
    plt.figure(figsize=(8,5))
    plt.plot(times, diffusion)
    plt.xlabel("Time")
    plt.ylabel("Diffusion")
    plt.title("Evolution of the diffusion")
    plt.grid()
    plt.show()

########################################################   
################# Specific Heat ########################
########################################################

def plot_specific_heat(total_specific_heat, particle_specific_heat, time_step):
    """
    This functions creates a plot of the diffusion.
    
    Parameters
    ----------    

    Total_specific_heat: np.array
        Total specific heat values over time
        Array dimension: number_of_time_steps
    Particle_specific_heat: np.array
        Total specific heat values over time
        Array dimension: number_of_time_steps
    time_step: float
        The time step used in the simulation
        
    Returns  
    ----------   
    """
    
    fig, ax = plt.subplots(1,2)
    number_time_steps=total_specific_heat.shape[0]
    times = np.arange(0, number_time_steps*time_step, time_step)

    ax[0].plot(times, total_specific_heat[~np.isnan(total_specific_heat)])
    
    # Plot vertical lines for NaN values
    nan_indices = np.where(np.isnan(total_specific_heat))[0]
    ax[0].vlines(nan_indices, ymin=min(total_specific_heat), ymax=max(total_specific_heat), color='r', linestyle='--')

    ax[1].plot(times, particle_specific_heat[~np.isnan(particle_specific_heat)])

    # Plot vertical lines for NaN values
    nan_indices = np.where(np.isnan(particle_specific_heat))[0]
    ax[1].vlines(nan_indices, ymin=min(particle_specific_heat), ymax=max(particle_specific_heat), color='r', linestyle='--')

    ax[0].set_xlabel("Time")
    ax[1].set_xlabel("Time")

    ax[0].set_ylabel(r"$C_V$")
    ax[1].set_ylabel(r'$c_V$')
    
    ax[0].set_title("Total Specific heat")
    ax[1].set_title("Specific heat per atom")
    
    plt.tight_layout()
    plt.show()



    
    