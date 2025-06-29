"""
In this file, we define methods to estimate the errors on the different observables in the simulation.
"""

from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

########################################################   
######### Autocorrelation analysis #####################
########################################################

def neg_exponential_fit(t, tau):
    return np.exp(-t/tau)

def neg_exponential_fit_w_offset(t, tau, offset, a):
    return offset-a*np.exp(-t/tau)


def autocorrelation_analysis(number_time_steps, obs, obs_name, start, plot=False):
    """
    This function determines the autocorrelation function and correlation time of a given observable.
    
    Parameters
    ----------
    number_time_steps: int
        The total number of simulation steps
    obs: np.ndarray
        Array with the value of a observable at different times
        length = number_time_steps
    obs_name: str
        Name of the observable (used for the title of the plot)
    start: int
        Index at which the autocorrelation function should start to be computed
    plot: bool
    
    Returns
    -------
    autocorrelation: np.ndarray
        Array that contains the autocorrelation function at the discretized times
        length = number_time_steps
    correlation_time: float
        Correlation time (time_index specifically) that gives the time distance at which the data values are uncorrelated
    obs_error: float
        Error of the observable using autocorrelation method
    """

    # Generate the time indices 
    time_indices = np.linspace(0, number_time_steps - start - 1, number_time_steps - start)
    # Convert them to integers
    time_indices = time_indices.astype(int)

    # Select the values of the observable from start to the end
    obs = obs[start:]

    # Recalculate the number of time steps removing the s first elements
    number_time_steps = number_time_steps - start

    # Calculation of the autocorrelation function at discretized timesteps, cf. Lecture Notes, Week 5 
    a = np.array([(number_time_steps-t)*np.sum(obs[:number_time_steps-t]*obs[t:]) for t in time_indices]) 
    b = np.array([np.sum(obs[:number_time_steps-t])*np.sum(obs[t:]) for t in time_indices])
    c = np.array([(number_time_steps-t)*np.sum(obs[:number_time_steps-t]**2) - np.sum(obs[:number_time_steps-t])**2 for t in time_indices])
    d = np.array([(number_time_steps-t)*np.sum(obs[t:]**2) - np.sum(obs[t:])**2 for t in time_indices])

    # Find the indices where both c and d are non-zero
    valid_indices = np.where((c != 0) & (d != 0))[0]
    
    # Compute the autocorrelation only for the valid indices
    autocorrelation = (a[valid_indices] - b[valid_indices]) / (np.sqrt(c[valid_indices]) * np.sqrt(d[valid_indices]))

    if plot:
        # Plot the autocorrelation function to identify the region to do the exponential fit
        plt.figure(figsize=(8,5))
        plt.plot(start+valid_indices, autocorrelation, label = 'Simulation data')  
        plt.xlabel('Time indices')
        plt.ylabel('Autocorrelation')
        plt.title(f'Autocorrelation of {obs_name}')

    # We observe from the graph that the autocorrelation behaves well until approximately the time_threshold^th time
    # We fit the autocorrelation to a negative exponential in this region to obtain the correlation time
    
    time_threshold = 100   # Note: Important to change for each observable
    autocorrelation_for_fit = autocorrelation[valid_indices < time_threshold]
    times_for_fit = valid_indices[valid_indices < time_threshold]

    fitting_parameters, fitting_parameters_cov = curve_fit(neg_exponential_fit, times_for_fit, autocorrelation_for_fit)
    correlation_time = fitting_parameters[0]
    if plot:
        plt.plot(start+times_for_fit, neg_exponential_fit(times_for_fit, correlation_time), label = f'Exponential fit (correlation time = {correlation_time})')
        plt.legend()
        plt.grid()
        plt.show()

    obs_error = np.sqrt((2*correlation_time/(number_time_steps))*(np.mean(obs**2)-np.mean(obs)**2))
    
    return autocorrelation, correlation_time, obs_error


def normal_autocorr(mu, sigma, tau, N):
    """
    Generates an autocorrelated sequence of Gaussian random numbers.
    
    Each of the random numbers in the sequence of length `N` is distributed
    according to a Gaussian with mean `mu` and standard deviation `sigma` (just
    as in `numpy.random.normal`, with `loc=mu` and `scale=sigma`). Subsequent
    random numbers are correlated such that the autocorrelation function
    is on average `exp(-n/tau)` where `n` is the distance between random
    numbers in the sequence.
    
    This function implements the algorithm described in
    https://www.cmu.edu/biolphys/deserno/pdf/corr_gaussian_random.pdf
    
    Parameters
    ----------
    mu: float
        mean of each Gaussian random number
    sigma: float
        standard deviation of each Gaussian random number
    tau: float
        autocorrelation time
    N: int
        number of desired random numbers
    
    Returns
    --------
    sequence: numpy array
        array of autocorrelated random numbers
    """
    
    f = np.exp(-1./tau)
    
    sequence = np.zeros(shape=(N,))
    
    sequence[0] = np.random.normal(0, 1)
    for i in range(1, N):
        sequence[i] = f * sequence[i-1] + np.sqrt(1 - f**2) * np.random.normal(0, 1)
    
    return mu + sigma * sequence


########################################################   
#################### Data Blocking #####################
########################################################

def data_blocking(number_time_steps, start, obs, obs_name, plot, initial_guess):
    """
    We compute the error of an observable using data blocking

    Parameters
    ----------
    number_time_steps: int
        The total number of simulation steps
    start: int
        Starting index of the data to analyse
    obs: np.ndarray
        Observable to analyse
    obs_name: str
        Name of the observable to analyse
    plot: boolean
        Parameter to determine whether or not a plot with the Standard Deviations in terms of the block size should be made

    Returns
    --------
    sigma: np.array
        Array with all the standard deviations for all the values of the block size b.
    estimated_error_fit: float
        Computed error for the specified observable using the convergence criteria with an exponential fit.
    """
    
    total_time_steps=number_time_steps-start
    b_index=np.array([i for i in range(1, int(total_time_steps/5)+1)]) # here we assume that total_time_steps is large enough to see a convergence of the data
    Blocks = total_time_steps / b_index

    sigma=np.zeros_like(b_index, dtype=float)
    
    for index, b in enumerate(b_index):
        a=(1/b)*np.array([np.sum(obs[start+i*b: start+(i+1)*b]) for i in range(0, int(Blocks[index]))])

        mean_a1=np.mean(a)      # We compute the mean of a
        mean_a2=np.mean(a*a)    # We compute the mean of a**2

        variance=(mean_a2-mean_a1*mean_a1)*(1/(int(Blocks[index])-1))
        sigma[index]=np.sqrt(variance)

    if plot:
        # Plot the autocorrelation function to identify the region to do the exponential fit
        plt.figure(figsize=(8,5))
        plt.plot(b_index, sigma, label = 'Simulation data')
        plt.xlabel('Blocks b indices')
        plt.ylabel('$\sigma(b)$')
        plt.title(f'Standard deviation in terms of the block size for the observable:{obs_name}')
        

    #Let us get the error with an exponential fitting 
    fitting_parameters, fitting_parameters_cov = curve_fit(neg_exponential_fit_w_offset, b_index, sigma, p0=initial_guess)
    correlation_time = fitting_parameters[0]
    estimated_error_fit = fitting_parameters[1]
    amplitude = fitting_parameters[2]
    
    if plot:
        plt.plot(b_index, neg_exponential_fit_w_offset(b_index, correlation_time, estimated_error_fit, amplitude), label = f'Exponential fit (correlation time = {correlation_time})')
        plt.legend()
        plt.grid()
        plt.show()
    
    return sigma, estimated_error_fit


########################################################   
############### Bootstrapping ##########################
########################################################

def bootstrapping_msd(iterations, time_series_diffusion):
    """
    This function implements the bootstrapping method to analyze an error in the diffusion observable,
    for a certain time step, given access to the time series of diffusion.
    
    Parameters
    ----------
    iterations: int
        The number of bootstrap iterations
    time_series_diffusion: np.array
        Time series data for observable calculation
        
    Returns
    ----------
    error: float
        Error estimate for diffusion observable
    """
    
    # Preprocess data by applying data blocking
    values_per_block = 250
    number_of_blocks = int(time_series_diffusion.shape[0]/values_per_block) # 250 values per block
    diffusion_ts_blocked = np.array([1/values_per_block*np.sum(time_series_diffusion[i*values_per_block:(i+1)*values_per_block]) for i in range(number_of_blocks)])   
    
    iterations = iterations
    diffusion_bootstrap = np.zeros(iterations)
    
    for i in range(iterations): 
        diffusion_ts_bootstrap = np.random.choice(diffusion_ts_blocked, diffusion_ts_blocked.shape[0], replace=True) # drawing position vectors with replacement by sampling indices
        diffusion_mean_bootstrap = np.mean(diffusion_ts_bootstrap) # determine observable
        diffusion_bootstrap[i] = diffusion_mean_bootstrap          # append diffusion term
    
    error = np.sqrt( (np.mean(diffusion_bootstrap**2) - np.mean(diffusion_bootstrap)**2 ))
    return error


def bootstrapping_specific_heat(iterations, e_kin_ts, num_atoms):
    """
    This function implements the bootstrapping method to analyze an error in the specific heat observable for a certain time step,
    given access to the time series of kinetic energies.
    
    Parameters
    ----------
    iterations: int
        The number of bootstrap iterations

    e_kin_ts: np.array
        Time series data for observable calculation
        
    num_atoms: int
        Number of atoms in simulation
        
    Returns
    ----------
    error: float
        Error estimate for specific heat observable
    """
    
    iterations = iterations
    spec_heat_bootstrap = np.zeros(iterations)
    
    # Preprocess data by applying data blocking
    values_per_block = 250
    number_of_blocks = int(e_kin_ts.shape[0]/values_per_block) # 250 values per block
    e_kin_ts_blocked = np.array([1/values_per_block*np.sum(e_kin_ts[i*values_per_block:(i+1)*values_per_block]) for i in range(number_of_blocks)]) 
    
    for i in range(iterations):  
        # Bootstrap sampling of blocked kinetic energies
        e_kin_ts_bootstrap = np.random.choice(e_kin_ts_blocked, e_kin_ts_blocked.shape[0], replace=True) 
        
        # Determine total specific heat capacity from blocked data for each sample
        delta_K = np.mean(e_kin_ts_bootstrap**2)-np.mean(e_kin_ts_bootstrap)**2
        K_sq = np.mean(e_kin_ts_bootstrap)**2

        cv_inv = 2/(3*num_atoms)-delta_K/K_sq
        cv = np.nan if cv_inv == 0 else 1/cv_inv
        
        spec_heat_bootstrap[i] = cv # append specific heat term
    
    error = np.sqrt( (np.nanmean(spec_heat_bootstrap**2) - np.nanmean(spec_heat_bootstrap)**2 ))
    
    return error
