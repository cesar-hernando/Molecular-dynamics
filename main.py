"""
authors: Arturo, Cesar, Tim

In this file, we retrieve the parameters and iteratively simulate the system, plotting the evolution and
storing observables along the way.

Credit: References include Joseph Thijssen's Computational Physics book and the Computational Physics Lecture Notes, accessible at https://compphys.quantumtinkerer.tudelft.nl.
"""

import parameters, simulation, animation, plots, units_conversion, error_estimation, observables
import numpy as np


###############################################################################
####################### Retrieving parameters #################################
###############################################################################

rho, num_atoms_prior, time_step, number_time_steps = parameters.rho, parameters.num_atoms, parameters.time_step, parameters.number_time_steps
temperature, rescaling, rescaling_duration = parameters.Normal_Temp, parameters.rescaling, parameters.rescaling_duration
method = parameters.method
bootstrap_iterations = parameters.bootstrap_iterations

box_dimension = (num_atoms_prior/rho)**(1/3)
scaling_factor = np.round((num_atoms_prior/4)**(1/3))  # determine factor for lattice constant scaling
lattice_constant = box_dimension/scaling_factor        # determine lattice constant from box dimension and scaling factor

dimensionality = parameters.dimensionality
bound=parameters.boundary_representation

###############################################################################
####################### Initialization ########################################
###############################################################################

init_positions, num_atoms = simulation.fcc_lattice(box_dimension, lattice_constant)  # Retrieve initial positions and velocities of the Argon atoms

assert num_atoms_prior == num_atoms

init_velocities = simulation.init_velocity(num_atoms, temperature)

print("\n")
print(f"Argon Molecular Simulation of {num_atoms} atoms in a box of dimension {box_dimension:.4f} x {box_dimension:.4f} x {box_dimension:.4f}.\n __________________________________ \n")

# Plot initial velocities
plots.velocity_distribution(init_velocities, temperature)                      



###############################################################################
####################### Simulation ############################################
###############################################################################

# Perform simulation
positions, velocities, kinetic_energy, potential_energy, distances, Temperatures = simulation.simulate(init_positions, init_velocities, number_time_steps, time_step, box_dimension, num_atoms, method, rescaling, temperature, rescaling_duration)

# Calculate the total energy
total_energy = kinetic_energy + potential_energy

# Analyse autocorrelation of kinetic energy and calculate error
start=int(number_time_steps*rescaling_duration)
autocorrelation_kin_en, tau_kinetic_energy, error_kinetic_energy = error_estimation.autocorrelation_analysis(number_time_steps, kinetic_energy, 'Kinetic energy', start=start, plot=True)

print("\n")
print(f"Total Kinetic Energy = {np.mean(kinetic_energy)} +- {error_kinetic_energy} (autocorrelation)")

# Analyse autocorrelation of potential energy and calculate error
autocorrelation_pot_en, tau_potential_energy, error_potential_energy = error_estimation.autocorrelation_analysis(number_time_steps, potential_energy, 'Potential energy', start=start, plot=True)
print(f"Total Potential Energy = {np.mean(potential_energy)} +- {error_potential_energy} (autocorrelation) \n")

# Plot energies
plots.plot_energies(number_time_steps, time_step, kinetic_energy, potential_energy, total_energy)

# Plot temperatures
plots.plot_temp(number_time_steps, time_step, Temperatures, temperature, SI=False)


###############################################################################
###################### Derived Observables ####################################
###############################################################################

                #### Determine diffusion and mean-squared displacement ####

mean_squared_displacement, diffusion = observables.diffusion(positions, box_dimension, time_step, number_time_steps = number_time_steps, initial_time_index = start)
plots.plot_diffusion(time_step, mean_squared_displacement, diffusion)                # Plot diffusion

#Obtain the error using Data Blocking
sigma_MSD, error_DFS_DB = error_estimation.data_blocking(number_time_steps, start=start, obs=diffusion, obs_name='Diffusion', plot=True, initial_guess=[10.0, 0.5, 0.5])

time_avg_diffusion=np.mean(diffusion)

print(f"Diffusion = {time_avg_diffusion} +- {error_DFS_DB} (Data blocking)\n")

# Error estimate for time-averaged mean square displacement
time_avg_mean_squared_displacement = np.mean(mean_squared_displacement)
msd_error = error_estimation.bootstrapping_msd(bootstrap_iterations, mean_squared_displacement[(int(number_time_steps*rescaling_duration)):])

#Obtain the error using Data Blocking
sigma_MSD, error_MSD_DB = error_estimation.data_blocking(number_time_steps, start=start, obs=mean_squared_displacement, obs_name='Mean squared displacement', plot=False, initial_guess=[10.0, 0.5, 0.5])

print(f"Time average mean squared displacement = {time_avg_mean_squared_displacement} +- {msd_error} (Bootstrap)")
print(f"Time average mean squared displacement = {time_avg_mean_squared_displacement} +- {error_MSD_DB} (Data blocking)\n")


                #### Determine specific heat ####

Total_specific_heat, Particle_specific_heat = observables.Specific_Heat(Kin_En=kinetic_energy, num_atoms = num_atoms, number_time_steps = number_time_steps, initial_time = int(number_time_steps * rescaling_duration) )
# The initial time step to start computing the specific heat is once the rescaling is finished, i.e. number_time_steps * rescaling_duration

# Plot specific heat
plots.plot_specific_heat(Total_specific_heat, Particle_specific_heat, time_step=time_step)

# Error estimate for specific heat with Bootstrap
heat_error_bootstrap = error_estimation.bootstrapping_specific_heat(bootstrap_iterations, kinetic_energy[(int(number_time_steps*rescaling_duration)):], num_atoms)
heat_sigmas_DataBlock, heat_error_DB = error_estimation.data_blocking(number_time_steps-int(number_time_steps * rescaling_duration), start=0, obs=Total_specific_heat, obs_name='Specific_Heat', plot=True, initial_guess=[10.0, 0.5, 0.5])

print(f"Specific Heat/Particle at end of simulation: {Particle_specific_heat[int(number_time_steps *(1-rescaling_duration))-1]}")
print(f"Total Specific Heat at end of simulation: {Total_specific_heat[int(number_time_steps *(1-rescaling_duration))-1]} +- {heat_error_bootstrap} (Bootstrap)")
print(f"Total Specific Heat  at end of simulation:  {Total_specific_heat[int(number_time_steps *(1-rescaling_duration))-1]} +- {heat_error_DB} (Data Blocking)\n")


                #### Pair correlation ####

bin_partition, correlation = observables.pair_correlation(distances[:,:,(int(number_time_steps*rescaling_duration)):], box_dimension, num_atoms)
plots.plot_pair_correlation(bin_partition, correlation)



###############################################################################
####################### Animation #############################################
###############################################################################
 
plot = animation.animate(positions, box_dimension, dimensionality, bound)


###############################################################################
##################### Observables in SI units #################################
###############################################################################

rho_SI=units_conversion.dimensionless_to_SI(rho, "number_density")
print(f"In SI units, the considered number density of the system is {rho_SI} m^(-3) ")
Temperatures_SI=units_conversion.dimensionless_to_SI(Temperatures,  "temperature")
Objective_T_SI=units_conversion.dimensionless_to_SI(temperature,  "temperature")
plots.plot_temp(number_time_steps, time_step, Temperatures_SI, Objective_T_SI, SI=True)

Total_Specific_Heat_SI_last=units_conversion.dimensionless_to_SI(Total_specific_heat[int(number_time_steps *(1-rescaling_duration))-1], "specific_heat")
Total_heat_error_SI=units_conversion.dimensionless_to_SI(max(heat_error_bootstrap, heat_error_DB), "specific_heat")

Particle_Specific_Heat_SI_last=units_conversion.dimensionless_to_SI(Particle_specific_heat[int(number_time_steps *(1-rescaling_duration))-1], "specific_heat")
Particle_heat_error_SI=Total_heat_error_SI/num_atoms
print(f"The value of the total specific heat in SI units: {Total_Specific_Heat_SI_last} +- {Total_heat_error_SI} J·K^(-1)·kg^(-1).")
print(f"The value of the specific heat per particle in SI units: {Particle_Specific_Heat_SI_last} +- {Particle_heat_error_SI} J·K^(-1)·kg^(-1).\n")


###############################################################################
####################### Testing ##############################################
###############################################################################

# Test implementation of autocorrelation function
tau_random_data = 50
number_data_points = 1000
random_gen_data = error_estimation.normal_autocorr(0, 1, tau_random_data, number_data_points) 
random_data_autocorr, estimated_tau, random_obs_error = error_estimation.autocorrelation_analysis(number_data_points, random_gen_data, 'random sequence', start=0, plot=True)
estimated_gen_error = np.sqrt((2*estimated_tau/(number_data_points))*(np.mean(random_gen_data**2)-np.mean(random_gen_data)**2))
print(f"Imposed tau = {tau_random_data}, Estimated tau = {estimated_tau:2f}, Estimated error: {estimated_gen_error:2f} (autocorrelation)")

# Test implementation of Data Blocking
number_data_points = 1000
random_gen_data = error_estimation.normal_autocorr(0, 1, tau_random_data, number_data_points) 
sigma, error_fit_DB=error_estimation.data_blocking(number_data_points, start=0, obs=random_gen_data, obs_name='random sequence', plot=True, initial_guess=[1.0, 0.5, 0.5])
print(f"Estimated error in Randomly generated data:{error_fit_DB} (DataBlocking)")



