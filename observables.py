"""
In this file the observables are computed.
"""

import numpy as np

##################################################################
###################### Energies ##################################
##################################################################

def kinetic_energy(vel):
    """
    Computes the dimensionless kinetic energy of an atomic system.

    Parameters
    ----------
    vel: np.ndarray of dimension 3 x N
        Velocity of particle

    Returns
    -------
    float
        The total kinetic energy (dimensionless) of the system.
    """
    
    E_kin = 0.5*np.sum(vel**2) 
     
    return E_kin


def potential_energy(rel_dist):
    """
    Computes the potential energy of an atomic system.

    Parameters
    ----------
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances
        Array dimension: N x N

    Returns
    -------
    float
        The total potential energy of the system.
    """

    # We retrieve the indices from the upper triangle of the matrix without considering the diagonal (k=0) in 
    # order to avoid double counting terms
    
    ind_up_triangle=np.triu_indices_from(rel_dist, k=1) 
    
    # Using these indices and having access to the N x N array of relative distances, we determine the potential energy 
    # of the system
    
    U=np.sum((rel_dist[ind_up_triangle])**(-12)-(rel_dist[ind_up_triangle])**(-6))
    U = 4 * U
    
    return U


##################################################################
################### Pair correlation #############################
##################################################################

def pair_correlation(distances, box_dim, num_atoms):
    """
    Computes the pair correlation function of the atomic system.

    Parameters
    ----------
    distances : np.ndarray
        Relative particle distances as obtained from atomic_distances
        Array dimension: N x N x number_time_steps (post equilibration)
    box_dim : float
        The dimension of the simulation box
    num_atoms : int
        The number of particles in the simulation

    Returns
    -------
    correlation: np.ndarray
        The pair correlation function (array containing function values)
    """

    # Define bins, i.e., a domain partition, for the pair correlation function
    num_bins = 250
    
    # Initialize correlation function and create a partition of the domain
    correlation = np.zeros(num_bins)
    bin_partition = np.linspace(0, np.sqrt(3/4*box_dim**2), num_bins+1) # Note that the maximal distance is upper bounded by \sqrt{3/4*box_dim^2} due to PBC
    delta_r = np.sqrt(3/4*box_dim**2) / num_bins 

    # Compute the histogram using numpy functionality
    distance_histogram, bin_partition = np.histogram(distances, bins=bin_partition)
    
    # Discard the self-distances and divide all others by 2 to account for overcounting pairs
    distance_histogram[0] = 0
    distance_histogram = distance_histogram / 2
    bin_partition = bin_partition[:-1]

    # Calculating pairwise correlation function (apply normalization factor) per bin using numpy vectorization
    correlation = 1/distances.shape[2] * (2*box_dim**3)/(num_atoms*(num_atoms-1)) * distance_histogram/(4*np.pi*(bin_partition+10**(-12))**2*delta_r ) 
  
    assert correlation.shape[0] == num_bins
    
    return bin_partition, correlation


    
##################################################################
########################## Diffusion #############################
##################################################################

def diffusion(positions, box_dim, timestep, number_time_steps=None, initial_time_index=None):
    """
    Computes the average (with respect to the particles) standard deviation of the positions with respect to an initial equilibrium position. 

    Parameters
    ----------
    positions: np.ndarray 
        Positions of the atoms at each time step
        array dimensions: num_atoms x 3 x number_time_steps
    box_dim : float
        The dimension of the simulation box
    timestep : float
        Duration of a single simulation step
    number_time_steps: int
        The total number of simulation steps
    initial_time_index: int
        Index of the initial time considered to use as a reference point in order to calculate the standard deviation of positions
        
    Returns
    -------
    mean_squared_displacement: np.ndarray
        Average square distance of particles with respect to their positions at initial_time.
        Formula: <|r(t)-r(0)|^2>
        array length: number_time_steps
    diffusion: np.ndarray
        Diffusion coefficient given by the formula: mean_squared_displacement/(6t)
        array length: number_time_steps
    """

    # Obtain number of time steps if it has not been indicated as an input
    if number_time_steps is None:
        number_time_steps = positions.shape[2]

    # Choose index of time to consider for initial positions if it has not been indicated as an input
    if initial_time_index is None:
        initial_time_index = number_time_steps//3

    # Calculate square distances between positions at different times and the initial time, taking into account perioduc boundary conditions

    # 1. Remove positions prior to initial time considered
    positions = positions[:,:,initial_time_index:]
    # 2. Select initial positions
    initial_positions = positions[:,:,0]
    # 3. Calculate the displacement of each particle with respect to their initial (equilibrium) position
    displacement = positions - initial_positions[:,:,np.newaxis]
    # 4. In order to account for the minimal image convention, we apply the following formula retrospectively.
    displacement -= (displacement + box_dim/2) // (box_dim) * box_dim
    # 5. Calcultae the norm squared of the displacement
    squared_displacement = np.sum(displacement*displacement, axis=1)
    # 6. Average over all the particles
    mean_squared_displacement = np.mean(squared_displacement, axis=0)

    # 7. Calculate the diffusion
    times = np.arange(timestep*initial_time_index, number_time_steps*timestep, timestep)
    diffusion = mean_squared_displacement/(6*times)

    return mean_squared_displacement, diffusion


##################################################################
###################### Specific heat #############################
##################################################################

def Specific_Heat(Kin_En, num_atoms, number_time_steps, initial_time):
    """
    Computes the specific heat (total and per particle) at each time step after the scaling has ended. 

    Parameters
    ----------
    Kin_En: np.ndarray 
        Kinetic energy of the system at each time step
        array dimensions: number_time_steps (1D)
    num_atoms : int
        Number of particles in the system
    number_time_steps: int
        The total number of simulation steps
    initial_time_index: int
        Index of the initial time considered to use as a reference point in order to calculate the standard deviation of positions
        
    Returns
    -------
    Total_specific_heat : np.ndarray
        Total specific heat at each time step greater than initial_time_steps.
        Formula: [1/(3N)+1/2-(<K**2>/(2<K>**2))]**(-1)
        array length: number_time_steps-initial_time_steps

    Particle_specific_heat : np.ndarray
        Specific heat per particle at each time step greater than initial_time_steps.
        Formula: (Total_specific_heat)/num_atoms
        array length: number_time_steps-initial_time_steps
    """

    #Let us first compute the averages of K and K^2
    Cumulation_Kin=np.cumsum(Kin_En[initial_time:])
    Average_Kin=Cumulation_Kin/(np.arange(1, number_time_steps-initial_time+1)) #<K>=(1/(n-n0))*sum(K_i, i=initial_time_steps+1=n0+1, n), where the last two elements are the beginning and end of the sum

    Cumulation_Kin2=np.cumsum(Kin_En[initial_time:]*Kin_En[initial_time:])
    Average_Kin2=Cumulation_Kin2/(np.arange(1, number_time_steps-initial_time+1)) #<K^2>=(1/(n-n0))*sum(K_i^2, i=initial_time_steps+1=n0+1, n)

    #We now compute <K>^2
    Square_Average_Kin=Average_Kin*Average_Kin

    Inverse_Cv=2/(3*num_atoms)+1-(Average_Kin2)/(Square_Average_Kin)
    mask= Inverse_Cv!=0

    Total_specific_heat=np.full_like(Inverse_Cv, np.nan, dtype=float) #We initialise it full of NaN. This way, if there is any value in Inverse_Cv that is equal to 0 we save it as a NaN.
    Total_specific_heat[mask]=1/Inverse_Cv[mask]

    Particle_specific_heat=Total_specific_heat/num_atoms

    return Total_specific_heat, Particle_specific_heat

    