"""
In this file, we plot the state of the system for a given time.
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

def animate(positions, box_dim, dimensionality = 3, boundary = False):
    """
    Molecular dynamics animation. Plots positions of the Argon atoms.

    Parameters
    ----------
    positions : np.ndarray
        The positions of the atoms to be plotted, in Cartesian space
    box_dim : float
        Dimensions of the simulation box
    dimensionality: int
        Flag for dimensionality of animated system
    boundary: bool
        Determines if boundary boxes are plotted or not
    
    Returns
    -------
    None
    """
    
    ########################################################   
    ############### 2D system ##############################
    ########################################################
    
    if dimensionality == 2:
        
        # Create a plot with nine boxes
        fig, ax = plt.subplots(figsize=(9,9))                       
        plt.axhline(-box_dim/2, linestyle='--', color='black')
        plt.axhline(box_dim/2, linestyle='--', color='black')
        plt.axvline(-box_dim/2, linestyle='--', color='black')
        plt.axvline(box_dim/2, linestyle='--', color='black')
        
        plt.xlim(-3*box_dim/2,3*box_dim/2)
        plt.ylim(-3*box_dim/2,3*box_dim/2)
        
        # Create empty scatter plot
        sc = ax.scatter([], [])                                     
        
        # Animation update function
        def update(frame):                                          
            x = positions[:, 0, frame]
            y = positions[:, 1, frame]

        # Create periodic copies of the particles
            copies_x = np.concatenate([                             
            x - box_dim, x, x + box_dim,
            x - box_dim, x, x + box_dim,
            x - box_dim, x, x + box_dim,
            ])
        
            copies_y = np.concatenate([
            y + box_dim, y + box_dim, y + box_dim,
            y, y, y,
            y - box_dim, y - box_dim, y - box_dim,
            ])

            all_positions = np.c_[copies_x, copies_y] 
            sc.set_offsets(all_positions)
            return sc,
        
        ani = animation.FuncAnimation(fig, update, frames=positions.shape[2], interval=5, blit=True)   # Create animation

        plt.show()

    ########################################################   
    ############### 3D system with boundary ################
    ########################################################
    
    elif dimensionality == 3 and boundary==True:

        fig = plt.figure(figsize=(9,9))
        ax = fig.add_subplot(projection='3d')
        
        ax.set_xlim(-3*box_dim/2, 3*box_dim/2)          
        ax.set_ylim(-3*box_dim/2, 3*box_dim/2)
        ax.set_zlim(-3*box_dim/2, 3*box_dim/2)
        
        # Shift copies along the xy-plane
        ax.plot([-box_dim/2,-box_dim/2], [-3*box_dim/2, 3*box_dim/2], [-box_dim/2, -box_dim/2], color='gray')  
        ax.plot([box_dim/2, box_dim/2], [-3*box_dim/2, 3*box_dim/2], [-box_dim/2, -box_dim/2], color='gray')  
        ax.plot([-3*box_dim/2, 3*box_dim/2], [box_dim/2, box_dim/2], [-box_dim/2, -box_dim/2], color='gray')  
        ax.plot([-3*box_dim/2, 3*box_dim/2], [-box_dim/2, -box_dim/2], [-box_dim/2, -box_dim/2], color='gray')  
        
        ax.plot([-box_dim/2,-box_dim/2], [-3*box_dim/2, 3*box_dim/2], [box_dim/2, box_dim/2], color='gray')  
        ax.plot([box_dim/2, box_dim/2], [-3*box_dim/2, 3*box_dim/2], [box_dim/2, box_dim/2], color='gray')  
        ax.plot([-3*box_dim/2, 3*box_dim/2], [box_dim/2, box_dim/2], [box_dim/2, box_dim/2], color='gray')  
        ax.plot([-3*box_dim/2, 3*box_dim/2], [-box_dim/2, -box_dim/2], [box_dim/2, box_dim/2], color='gray') 
        
        # Shift copies along the xz-plane
        ax.plot([-box_dim/2,-box_dim/2], [-box_dim/2, -box_dim/2], [-3*box_dim/2, 3*box_dim/2], color='gray')  
        ax.plot([box_dim/2, box_dim/2], [-box_dim/2, -box_dim/2], [-3*box_dim/2, 3*box_dim/2], color='gray')  
        ax.plot([-3*box_dim/2, 3*box_dim/2], [-box_dim/2, -box_dim/2], [box_dim/2, box_dim/2], color='gray')  
        ax.plot([-3*box_dim/2, 3*box_dim/2], [-box_dim/2, -box_dim/2], [-box_dim/2, -box_dim/2], color='gray')  
        
        ax.plot([-box_dim/2,-box_dim/2], [box_dim/2, box_dim/2], [-3*box_dim/2, 3*box_dim/2], color='gray')  
        ax.plot([box_dim/2, box_dim/2], [box_dim/2, box_dim/2], [-3*box_dim/2, 3*box_dim/2], color='gray')  
        ax.plot([-3*box_dim/2, 3*box_dim/2], [box_dim/2, box_dim/2], [box_dim/2, box_dim/2], color='gray')  
        ax.plot([-3*box_dim/2, 3*box_dim/2], [box_dim/2, box_dim/2],[-box_dim/2, -box_dim/2],  color='gray')   

        # Create empty 3D scatter plot
        sc = ax.scatter([], [], [])  

        # Define animation update function
        def update_3d(frame):
            x = positions[:, 0, frame]
            y = positions[:, 1, frame]
            z = positions[:, 2, frame]
            
        # Create periodic copies of the particles for 3D
            copies_x = np.concatenate([
                x - box_dim, x, x + box_dim,
                x - box_dim, x, x + box_dim,
                x - box_dim, x, x + box_dim,
                
                x - box_dim, x, x + box_dim,
                x - box_dim, x, x + box_dim,
                x - box_dim, x, x + box_dim,
                
                x - box_dim, x, x + box_dim,
                x - box_dim, x, x + box_dim,
                x - box_dim, x, x + box_dim,
            ])
        
            copies_y = np.concatenate([
                y + box_dim, y + box_dim, y + box_dim,
                y, y, y,
                y - box_dim, y - box_dim, y - box_dim,
                
                y + box_dim, y + box_dim, y + box_dim,
                y, y, y,
                y - box_dim, y - box_dim, y - box_dim,
                
                y + box_dim, y + box_dim, y + box_dim,
                y, y, y,
                y - box_dim, y - box_dim, y - box_dim,
            ])
            
            copies_z = np.concatenate([
                z - box_dim, z - box_dim, z - box_dim,
                z - box_dim, z - box_dim, z - box_dim,
                z - box_dim, z - box_dim, z - box_dim,
                
                z, z, z,
                z, z, z,
                z, z, z,

                z + box_dim, z + box_dim, z + box_dim,
                z + box_dim, z + box_dim, z + box_dim,
                z + box_dim, z + box_dim, z + box_dim,
            ])

            all_positions = np.c_[copies_x, copies_y, copies_z]
            sc._offsets3d = (all_positions[:, 0], all_positions[:, 1], all_positions[:, 2])
            return sc,

        # Call animation
        ani = animation.FuncAnimation(fig, update_3d, frames=positions.shape[2], interval=5, blit=False)

        plt.show()

    ########################################################   
    ############ 3D system without boundary ################
    ########################################################
    
    elif dimensionality == 3 and boundary==False:

        fig = plt.figure(figsize=(9,9))
        ax = fig.add_subplot(projection='3d')
        
        # Set limits for the plot
        ax.set_xlim(-box_dim/2, box_dim/2)
        ax.set_ylim(-box_dim/2, box_dim/2)
        ax.set_zlim(-box_dim/2, box_dim/2)
        
        ax.set_xticks(np.arange(-box_dim/2, box_dim/2, 1)) 
        ax.set_yticks(np.arange(-box_dim/2, box_dim/2, 1))
        ax.set_zticks(np.arange(-box_dim/2, box_dim/2, 1))

        # Create a scatter plot
        sc = ax.scatter([], [], [])  # Create empty 3D scatter plot

        # Animation update function
        def update_3d(frame):
            x = positions[:, 0, frame]
            y = positions[:, 1, frame]
            z = positions[:, 2, frame]

            all_positions = np.c_[x, y, z]
            sc._offsets3d = (all_positions[:, 0], all_positions[:, 1], all_positions[:, 2])
            return sc,
        
        # Create animation
        ani = animation.FuncAnimation(fig, update_3d, frames=positions.shape[2], interval=4, blit=False)

        plt.show()
        
    else:
        print("ERROR: Dimensionality not supported.")
        
