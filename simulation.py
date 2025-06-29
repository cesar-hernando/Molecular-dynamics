"""
In this file, we compute the time evolution of the system and the observables of interest.
"""

import numpy as np
import observables

#######################################################
############## MAIN SIMULATION ########################
#######################################################

def simulate(positions_init, velocities_init, number_time_steps, time_step, box_dim, num_atoms, method, rescaling, temp, rescaling_duration):
    """
    Molecular dynamics simulation using Euler or Verlet's algorithms
    to integrate the equations of motion. Calculates positions, velocities,
    energies and other observables at each timestep.
    
    Parameters
    ----------
    positions_init : np.ndarray
        The initial positions of the atoms in Cartesian space
        array dimensions: 3 x num_atoms
    velocities_init : np.ndarray
        The initial velocities of the atoms in Cartesian space
        array dimensions: 3 x num_atoms
    number_time_steps : int
        The total number of simulation steps
    time_step : float
        Duration of a single simulation step
    box_dim : float
        Dimensions of the simulation box
    num_atoms: float
        Number of atoms in the system
    method : str
        Can be either 'Euler', 'Verlet', 'RK4'
    rescaling: boolean
        Flag to indicate whether or not to use rescaling
    temp: float
        Target temperature for rescaling
    rescaling_duration: float
        Duration for which rescaling of system is applied 

    Returns
    -------
    positions: np.ndarray 
        Positions of the atoms at each time step
        array dimensions: num_atoms x 3 x number_time_steps
    velocities: np.ndarray 
        Velocities of the atoms at each time step
        array dimensions: num_atoms x 3 x number_time_steps
    distances: np.ndarray
        Relative distances between the atoms at every time step
        array dimensions: num_atoms x num_atoms x number_time_steps
    E_kin: np.ndarray 
        Kinetic energy at each time step
        array dimensions: 1 x number_time_steps
    E_pot: np.ndarray 
        Kinetic energy at each time step
        array dimensions: 1 x number_time_steps 
    """

    # Calculate the full simulation matrix, starting with some initial positions and velocities
    positions = np.zeros((num_atoms, 3, number_time_steps))         # The structure of the tensors is: (num_atoms, num_dimensions, num_time_steps), such that we have num_atoms layers of num_dimensions x num_time_steps matrices
    positions[:,:,0] = positions_init                               # set initial positions

    distances = np.zeros((num_atoms, num_atoms, number_time_steps)) # Here we save the distances between atoms 
    
    velocities = np.zeros((num_atoms, 3, number_time_steps))
    velocities[:,:,0] = velocities_init                             # set initial velocities
    
    E_kin = np.zeros((number_time_steps))
    E_pot = np.zeros((number_time_steps))

    Temp=np.zeros((number_time_steps))
    
    # Initial kinetic and potential energies
    E_kin[0] = observables.kinetic_energy(velocities[:,:,0])
    E_pot[0] = observables.potential_energy(atomic_distances(positions[:,:,0], box_dim=box_dim)[1])

    Temp[0]=temperature(num_atoms=num_atoms, kinetic=E_kin[0])

    # We simulate the evolution of the system with Euler's or Velocity-Verlet's method

    # Using Euler's method
    if method == "Euler":   
    
        for t in range(number_time_steps-1):
            # Retrieve positions and distances 
            rel_pos, rel_dist = atomic_distances(positions[:,:,t], box_dim=box_dim)
            # Store the distances at each time step
            distances[:, :, t]= rel_dist                          
            
            # Update positions
            positions[:,:,t+1] = positions[:,:,t] + velocities[:,:,t] * time_step
            
            # Implement Boundary conditions
            positions[:,:,t+1] = np.where(positions[:,:,t+1] > box_dim/2, positions[:,:,t+1]-box_dim, positions[:,:,t+1])
            positions[:,:,t+1] = np.where(positions[:,:,t+1] < -box_dim/2, positions[:,:,t+1]+box_dim, positions[:,:,t+1])
            
            # Update velocities
            velocities[:,:,t+1] = velocities[:,:,t] + lj_force(rel_pos=rel_pos, rel_dist=rel_dist)
                 
            # Determine observables
            E_kin[t+1] = observables.kinetic_energy(velocities[:,:,t+1])
            E_pot[t+1] = observables.potential_energy(rel_dist)

    
    # Using velocity Verlet's method:
    elif method == "Verlet":
        
        forces0=np.zeros(shape=(num_atoms, 3))
        forces1=np.zeros(shape=(num_atoms, 3))
        rescaling_count = 0
        
        for t in range(number_time_steps-1):
            if t % 250 == 0:
                print(f"Simulating time step {t}/{number_time_steps}: Simulation {np.floor(t/number_time_steps*100)} % done.")
            if t==0:
                # Retrieve positions and distances 
                rel_pos, rel_dist = atomic_distances(positions[:,:,t], box_dim=box_dim)
                distances[:, :, t]= rel_dist #Store the initial distances
    
                #Get F(x(t))
                forces0=lj_force(rel_pos=rel_pos, rel_dist=rel_dist)
                
            else:
                forces0=forces1 #We store the previous F(x(t+h)) as our new F(x(t))
    
            # Update positions
            positions[:,:,t+1] = positions[:,:,t] + velocities[:,:,t] * time_step + ((time_step**2)/2)*forces0
            
            # Implement Boundary conditions
            positions[:,:,t+1] = np.where(positions[:,:,t+1] > box_dim/2, positions[:,:,t+1]-box_dim, positions[:,:,t+1])
            positions[:,:,t+1] = np.where(positions[:,:,t+1] < -box_dim/2, positions[:,:,t+1]+box_dim, positions[:,:,t+1])
            
            # Compute F(x(t+h)) and the new positions
            if t<number_time_steps-1:
                rel_pos,rel_dist=atomic_distances(positions[:,:,t+1], box_dim=box_dim)
                distances[:, :, t+1]= rel_dist #Store the new distances
                forces1=lj_force(rel_pos=rel_pos, rel_dist=rel_dist)
                
                # Update velocities
                velocities[:,:,t+1] = velocities[:,:,t] + (time_step/2)*(forces1+forces0)

                # Determine observables
                E_kin[t+1] = observables.kinetic_energy(velocities[:,:,t+1])
                E_pot[t+1] = observables.potential_energy(rel_dist)
                Temp[t+1]=temperature(num_atoms=num_atoms, kinetic=E_kin[t+1])
                
                if rescaling == True:
                    if t<number_time_steps * rescaling_duration: # We rescale every rescaling_period time-steps during the fraction of the simulation given by rescaling_duration
                        if Temp[t+1]>(1.05)*temp or Temp[t+1]<(0.95)*temp:
                            Lambda=scaling(num_atoms, temp, E_kin[t])
                            velocities[:, :, t+1]=Lambda*velocities[:, :, t+1]
                            rescaling_count += 1  
        print(f"{rescaling_count} Rescalings applied")
        
    return positions, velocities, E_kin, E_pot, distances, Temp 


###############################################################################
########################## INITIALIZATION #####################################
###############################################################################

def fcc_lattice(box_dim, lat_constant):
    """
    Initializes a system of atoms on a Bravais-fcc lattice. This is done
    since it defines the ground state configuration of argon.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system
    lattice_const : float
        The lattice constant for an fcc lattice

    Returns
    -------
    pos_vec : np.ndarray
        Array of particle coordinates
        Array dimension: N x 3 (N = num_atoms) 
    num_atoms: int
        Number of particles that are initialized
    """

    ratio = int(box_dim/lat_constant)

    # Define the primitive basis vectors
    
    a_1 = 0.5*lat_constant*np.array([1, 1, 0]) 
    a_2 = 0.5*lat_constant*np.array([0, 1, 1]) 
    a_3 = 0.5*lat_constant*np.array([1, 0, 1]) 

    # Initialise list that stores positions of atoms in the FCC lattice
    positions = []

    # Initialize the atoms placed in the FCC lattice
    for n_1 in range(-ratio, 2*ratio + 1):
        for n_2 in range(-n_1, 2*ratio - n_1 + 1):
            for n_3 in range(-min(n_1,n_2), 2*ratio - max(n_1, n_2) + 1):
                positions.append(n_1*a_1 + n_2*a_2 + n_3*a_3)

    shift=0.5*box_dim*np.array([1, 1, 1])
    positions+=-shift
    
    # Convert positions list into an np array
    positions = np.array(positions) 

    positions = positions[np.all(positions < box_dim/2, axis=1)] # discards all rows (particles) that are duplicate through the lense of periodic boundary conditions
    
    # Count the number of atoms placed in the FCC lattice
    num_atoms = positions.shape[0]

    return positions, num_atoms


def init_velocity(num_atoms, temp):
    """
    Initializes the system with Gaussian distributed velocities.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system.
    temp : float
        The (unitless) temperature of the system.

    Returns
    -------
    vel_vec : np.ndarray
        Array of particle velocities 
        array dimensions: num_atoms * 3
    """
    
    # Initializing velocities of particles according to normalized Gaussian distribution
    velocities = np.random.normal(loc=0, scale=np.sqrt(temp), size = (num_atoms, 3)) # We incorporate the normalised temperature in the Maxwell-Boltzmann model
    
    # Center the sampled velocities    
    velocities = velocities-velocities.mean(axis=0)
    
    return velocities


#######################################################
################# FORCES ##############################
#######################################################

def atomic_distances(pos, box_dim):
    """
    This function determines both the relative positions and the relative distances between atoms
    from the array of atom positions in 3D space.
    
    TL;DR: This function uses multi-dimensional numpy vectorization.
    
    We have access to the N x 3 array of particle positions in Cartesian space at a given timestep t.
    First, we determine the relative positions. In order to do so efficiently, we use the numpy feature of broadcasting. 
    Here, we embed the N x 3 array with a new dimension, which casts it as N x 3 x 1 tensor. 
    Additionally, we transpose the original position array and embed it with a new dimension to create a
    1 x 3 x N tensor. Performing a pointwise operation between those tensor broadcasts the result to a 
    N x 3 x N tensor. In particular, we subtract the second tensor from the first. In effect, this results in 
    determining the relative position between every triple of positions (as stored in the second dimension of 
    the first tensor) with all other triples of positions (as stored in the third dimension of the second tensor). 
    This works since while traversing along the first dimension of the first tensor, for each step, we traverse through the third
    dimension of the second tensor and perform the subtractions.
    The broadcasted N x 3 x N tensor then contains the relative positions of all N particles, in all three Cartesian components. 
    
    Lastly, we determine the relative distances between particles. Recall that we have access to the N x 3 x N tensor of relative
    positions. Therefore, we only need to determine the norm of the entries along the second dimension (axis 0) in order to obtain a reduced N x N
    array which contains the pairwise atomic distances.
    
    Parameters
    ----------
    pos : np.ndarray
        The positions of the particles in cartesian space at a certain timestep t
        Array dimension: N x 3 (N = num_atoms) 
    box_dim : float
        The dimension of the simulation box

    Returns
    -------
    rel_pos : np.ndarray
        Relative positions of particles
        Array dimension: N x 3 x N
    rel_dist : np.ndarray
        The distance between particles
        Array dimension: N x N
    """

    rel_pos = pos[:,:,np.newaxis]-pos.T[np.newaxis,:,:] 
    
    # In order to account for the minimal image convention, we apply the following formula retrospectively.
    rel_pos -= (rel_pos + box_dim/2) // (box_dim) * box_dim 
    
    rel_dist = np.linalg.norm(rel_pos, axis=1)

    return rel_pos, rel_dist


def lj_force(rel_pos, rel_dist):
    """
    This function calculates the net forces on each atom, as dictated by the Lennard-Jones force.

    TL;DR: This function uses multi-dimensional numpy vectorization.
    
    (1)
    Here, we want to determine the LJ force that acts on each particle, as induced by all other particles in the system.
    We have access to the N x 3 x N tensor of relative positions, and the N x N tensor of relative distances. Note that the
    force depends on both the relative positions and distances. To access them in a homogeneous way, we add a new dimension 
    to the array of relative distances, which shapes it as a N x 1 x N tensor (the difference in the second dimension will be
    handled via broadcasting).
    
    (2)
    Next, note that the diagonal terms of the distance matrix are zero, since they refer to the self-distances. We need to account
    for these values in the LJ calculation below, since relative distances appear in the denominator there. Hence, we create a mask
    that stores the locations of the zero entries. In order to allow for some threshold (for internal Python floating point handling),
    we mark all elements that are less than 10**(-12) with a mask (another N x 1 x N tensor).
    
    (3)
    Having stored those locations, we proceed to overwrite the zero entries with a 1 such that we can use some numpy vectorization.
    The next line changes the existing N x 1 x N tensor to hold all relative distances, where zero self distances are overwritten with a 1
    
    (4)
    Now we come to the core part of the LJ force calculation. Having access to the N x 1 x N tensor of relative distances and the 
    N x 3 x N tensor of relative positions, we use the componentwise multiplication provided by numpy to leverage vectorization. The below
    force_terms tensor is the resulting (broadcasted) N x 3 x N tensor resulting from this multiplication.

    (5)
    Recalling that we had overwritten some elements, we need to account for the zero self-distances elements, whose positions are stored
    in mask. This can easily be done with numpy by copying the N x 1 x N mask to all rows alongg the seocnd dimension of the force_terms
    tensor, and replacing all terms with a zero. This is the appropriate correction since the force terms are summed over all other particles,
    and not over the particles themselves (i.e., a particle does not exert any force on itself.)
    
    (6)
    Having prepared the desired N x 3 x N tensor that contains the force terms acting on each particle as induced by all other particles (in
    each of the three Cartesian coordinates), we now sum the contributions. This is achieved by using the numpy summation along the third dimension
    (2nd axis). The result of this operation is a N x 3 array of forces that encapsulates the forces that each atom experiences in each component.
    
    Parameters
    ----------
    rel_pos : np.ndarray
        Relative particle positions as obtained from atomic_distances
        Array dimension: N x 3 x N
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances
        Array dimension: N x N

    Returns
    -------
    np.ndarray
        The net forces acting on each particle (in each component) due to LJ potential of all other particles
        Array dimension: N x 3
    """
    
    rel_dist = rel_dist.T[:,np.newaxis,:]                                                       # (1)
    
    mask = np.isclose(rel_dist, 0, atol = 10**(-12))                                            # (2)
    
    rel_dist1_without_zero = np.where(mask, 1, rel_dist)                                        # (3)
    
    force_terms = 24*(2*rel_dist1_without_zero**(-14)-rel_dist1_without_zero**(-8)) * rel_pos   # (4)
    
    force_terms = np.where(np.repeat(mask, 3, axis=1), 0, force_terms)                          # (5)
    
    forces = np.sum(force_terms, axis = 2)                                                      # (6)
    
    return forces


###############################################################################
########################## RESCALING ##########################################
###############################################################################


def temperature(num_atoms, kinetic):
    """
    This function returns the temperature as an observable of the system.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system.
    kinetic : float
        The kinetic energy of the system

    Returns
    -------
    T : float
        Temperature of the system
    """
    
    T = (2/3) * (kinetic / (num_atoms-1))
    
    return T

def scaling(num_atoms, temp, E_kin):
    """
    This function rescales the system to a given temperature.
    
    Parameters
    ----------
    num_atoms : int
        The number of particles in the system.
    temp : float
        The desired temperature of the system
    velocities: np.ndarray
        Array of particle velocities 
        array dimensions: num_atoms * 3

    Returns
    -------
    Lambda : float
        Rescaling factor that is to be applied to the velocities.
    """
    
    #norm_vel = np.linalg.norm(velocities, axis=1)
    #V = np.sum(norm_vel**2)

    V= 2*E_kin
    Lambda = np.sqrt(3 * (num_atoms-1) * (temp/V))
    
    return Lambda

    
    

