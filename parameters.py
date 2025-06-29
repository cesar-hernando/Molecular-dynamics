""""
In this file, various parameters of the system, simulation and animation are defined.
"""

########################################################   
############### System parameters ######################
########################################################

num_atoms = 108               # number of particles in the system, should be chosen from {4M^3, M \in \mathbb{N} } = {..., 32, 108, 256, 500, 864, ...}
rho = 1.06                     # specifies the density of the atoms in the system
Normal_Temp = 0.827               # normalised temperature, which corresponds to k*T/epsilon, k the Boltzmann constant


########################################################   
############### Simulation parameters ##################
########################################################

time_step = 0.004             # time step of simulation, chosen according to [Thijsen]
number_time_steps = 5000      # number of time steps of simulation
rescaling = True              # choose application of rescaling (if set to True, applies rescaling if temperature is outside a threshold of desired temperature)
rescaling_duration = 0.1      # fraction of simulation for which rescaling is applied (should be between 0.05 and 0.3)
method = 'Verlet'             # method to determine new positions and velocities from previous ones, default 'Verlet' (can be either 'Verlet' or 'Euler')
bootstrap_iterations = 10000  # defines the number of iterations used for bootstrap error estimation


########################################################   
############### Animation parameters ###################
########################################################

dimensionality = 3              # dimensionality of the system, can be either 2 or 3
boundary_representation = False # determine plotting of auxiliary PBC boxes, default is False

########################################################   
############### Universal constants ####################
########################################################
 
kb = 1.38E-23                   # [kg/K]

### Argon parameters
mass = 6.6E-26                  # [kg]
sigma = 3.405E-10               # [m]
temperature = 119.8             # [K]
eps = temperature*kb            # [J]
