# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 13:12:34 2014

kgsimulate simulates a solution for the 2nd order Klein-Gordon wave equation 
using a total spectral method. 

The simulation is ran in k-space (fourier transform of real space), in this 
form the (spectral) solution completely diaganolizes and we get Nk^3 
uncoupled equations for the fourier coefficients of our Electric Field.

Because of this, time-stepping is quite trivial, and the calculation is very
quick (but speed falls off in third order, which is usual of 3d problems)

Current problems to be addressed:
    Boundaries! With a discrete # of k values our problem is periodic, i.e. the
    simulation is running in a "box" and the waves will reflect on the edges.

@author: Stephen Carr, Brian Dawes
"""

import numpy as np
from numpy import pi

"""
Here we define our global variables
"""

#timestep h
h = 0.05

#total amount of timesteps
tsteps = 20

#mesh spacing in k-space
dk = 1

#number of mesh points in k-space
#note, odd number means we have k = 0
Nk = 101

#find max/min k values as a global variable
k_max = dk*((Nk-1)/2)

#mass
m = 1

"""
Next two methods are for our source, f
"""

#fourier coefficients of our source in the d'th component (d = 3 is z, etc)
def f_k(t,kx,ky,kz,d):
    
    if d == 1:
        return -4*pi*kx
        
    if d == 2:
        return -4*pi*ky
        
    if d == 3:  
        return -4*pi*kz

#matrix of the fourier coefficients in k space
def f(t,d):
    
    f_arr = np.zeros((Nk,Nk,Nk))
    
    for kx in range(-k_max,k_max,dk):
        for ky in range(-k_max,k_max,dk):
            for kz in range(-k_max,k_max,dk):
                f_arr[kx,ky,kz] = f_k(t,kx,ky,kz,d)
                
    return f_arr


"""
spectral(d) runs our spectral simulation for the d'th componenet of the E field
"""

def spectral(d):
    
    """ 
    define variables for our time stepping process
    for now assume inital conditions are zero
    Also set up k and mass arrays for the time step
    """
    
    #initial conditions (phi_last should always = phi_curr I think)
    phi_curr = np.zeros((Nk,Nk,Nk))
    phi_last = phi_curr
    
    #k^2 term as an array
    k_arr = np.zeros((Nk,Nk,Nk))
    for kx in range(-k_max,k_max,dk):
        for ky in range(-k_max,k_max,dk):
            for kz in range(-k_max,k_max,dk):
                k_arr[kx,ky,kz] = kx**2 + ky**2 + kz**2
    
    solution = np.zeros((tsteps+1,Nk,Nk,Nk))
    
    #running the simulation, only 5 lines!
    for t in range(0,tsteps):            
        phi_next = 2*phi_curr - phi_last + (h**2)*(f(t*h,d) - (k_arr+m**2)*phi_curr)
        solution[t,:,:,:] = phi_next
        phi_last = phi_curr
        phi_curr = phi_next
        print "Time Step " + str(t) + " done."
    
    return solution
    
def main():
    
    #Run the simulation in the d = 1,2,3 cases (corresponds to component x,y,z)
    K_x = spectral(1)
    K_y = spectral(2)
    K_z = spectral(3)
    
    #save our results
    np.save("K_x", K_x)
    np.save("K_y",K_y)
    np.save("K_z",K_z)
    
    #save the parameters (makes plotting easier)
    params = np.array([h, tsteps, dk, Nk, k_max, m])
    np.save("params",params)
    
if __name__ == "__main__":
    main()