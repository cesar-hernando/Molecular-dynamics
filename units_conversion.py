"""
In this file, we convert physical magnitudes from dimensionless units to SI units.
"""
import parameters
import numpy as np
from numpy import sqrt

def dimensionless_to_SI(value, magnitude_type):
    """
    This function accepts as input the value of a magnitude (e.g. 2.3) and its type (e.g. 'energy')
    and returns the value of the magnitude in SI units
    
    Parameters
    ----------
    value: float
        Input value to be transformed
    magnitude_type: str
        Allowed magnitude_type:
        - 'distance'
        - 'time'
        - 'velocity'
        - 'energy'
        - 'temperature'
        - 'specific_heat'
        - 'number_density'
        
    Returns
    ----------
    """
    
    if magnitude_type.lower() == 'distance':
        return value*parameters.sigma
        
    elif magnitude_type.lower() == 'time':
        return value*np.sqrt(parameters.mass*parameters.sigma**2/parameters.eps)
        
    elif magnitude_type.lower() == 'velocity':
        return value*sqrt(parameters.eps/parameters.mass)
        
    elif magnitude_type.lower() == 'energy':
        return value*parameters.eps

    elif magnitude_type.lower() == 'temperature':
        return value*parameters.eps/parameters.kb

    elif magnitude_type.lower() == "specific_heat":
        mass=parameters.mass
        kB=parameters.kb
        return value*(kB/mass) #We return the value of the Specific Heat in SI units: J·K^(-1)·kg^(-1)

    elif magnitude_type.lower() == "number_density":
        s=parameters.sigma
        return value/(s*s*s)

    else:
        raise ValueError(f"Invalid magnitude_type '{magnitude_type}'. Allowed types: 'distance', 'time', 'velocity', 'energy', 'temperature'.")

    


